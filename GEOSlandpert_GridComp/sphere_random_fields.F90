!
! sphere_random_fields.F90
!
! Module for generating isotropic Gaussian random fields directly on the
! sphere using the SMERFS/Markov construction (Creasey & Lang 2018).
! Candidate initialization replacement: Riccati logarithmic derivatives.
! The benchmark's MT/Ziggurat path still depends on the supplied SMERFS library.
!
! Parallel strategy (m-mode distribution):
!
!   The n_m Fourier m-modes are split evenly across the ranks of a
!   node-local communicator, established once in build_sphere_filter and
!   stored in sphere_filter_type for reuse.  The same strip assignment is
!   used for both the filter construction and the Kalman walk.
!
!   build_sphere_filter (comm optional):
!     Each rank builds the covariance / state-space chain only for its own
!     m-mode strip j_lo:j_hi and writes the result directly into shared
!     MPI windows holding all_innov / all_trans.  An MPI_Win_fence then
!     makes the complete filter visible to every rank.  No broadcast of the
!     filter tables is needed. Covariance integration is distributed within
!     each node (the full construction is repeated on other nodes).
!
!   draw_sphere_realisation:
!     Each rank walks the Kalman filter for its strip j_lo:j_hi of m-modes
!     independently (no communication during the walk) and writes the
!     partial half-spectrum into a shared window res_cplx(n_m, nz).  A
!     fence ensures all ranks see the complete spectrum before the IRFFT.
!     The IRFFT output rows are also distributed: each rank computes
!     rfield_out(:, row_lo:row_hi) into the shared rfield window.  After a
!     final fence every rank returns with a complete rfield_out(nphi, nz).
!
!   When no comm is supplied both subroutines fall back to a serial path
!   using private allocatables and a single strip covering everything.
!
! Public interface:
!   sphere_random_fields    - abstract_random_fields implementation
!   sphere_random_fields_id - cache key for the land_pert generator map
!
! Memory layout note (inherited from SMERFS example1_f.F90):
!   Arrays passed to C kernels use reversed Fortran dimension order so
!   column-major Fortran byte layout matches C row-major layout.
!
! wjiang, 2026-07
!
! -----------------------------------------------------------------------

module sphere_random_fields_mod
  use abstract_random_fieldsMod, only: abstract_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  use MAPL_MathConstantsMod, only: PI => MAPL_PI_R8, MAPL_DEGREES_TO_RADIANS_R8

  use, intrinsic :: iso_c_binding, only : c_double, c_float, c_int32_t, c_int64_t, &
                                          c_double_complex, c_loc,       &
                                          c_f_pointer, c_ptr,            &
                                           c_sizeof
  use mpi
  use nr_ran2_gasdev, only: NRANDSEED, nr_ran2_2d, nr_gasdev
  use smerfs_interface, only: zigg_f
  use sphere_riccati_init_mod, only: riccati_mode_filter, calibrate_spherical_c2, spherical_mode_cutoff
#ifdef MKL_AVAILABLE
  use MKL_DFTI
#endif

  implicit none
  private

  public :: sphere_random_fields
  public :: sphere_random_fields_id

  ! Filter order for the quadratic-in-Laplacian inverse spectrum.
  ! State: [g_m, (1-z*z)*dg_m/dz/kappa_m]. The Riccati initializer
  ! consistently scales the derivative; latitude reflection still flips its sign.
  integer, parameter :: MORD = 2
  ! Upper bound on omitted continuous-spectrum point variance (0.1%).
  ! This controls positive-m cutoff, independently of the output FFT size.
  real(c_double), parameter :: MODE_TAIL_TOL = 1.0d-3

  ! -----------------------------------------------------------------------
  !
  ! sphere_filter_type
  !
  ! Holds all pre-computed state-space filter coefficients for one
  ! choice of (nz, nphi, c0, c2), plus the MPI context needed by
  ! draw_sphere_realisation so it does not need to be re-established
  ! on every call.

  type :: sphere_filter_type

    ! Grid dimensions
    integer :: nz       = 0   ! total iso-latitude rings  (S.pole to N.pole)
    integer :: nphi     = 0   ! phi points per ring
    integer :: n_m      = 0   ! positive Fourier modes, including aliases
    integer :: uhalf    = 0   ! nz/2 + 1    (upper hemisphere ring count)

    ! Power-spectrum coefficients:  C_l = 1 / (c0 + c2*(l*(l+1))^2)
    real(c_double) :: c0 = 1.0d0
    real(c_double) :: c2 = 0.0d0

    ! Pre-computed state-space filter arrays (all m-modes):
    !   all_innov(MORD, MORD, uhalf,   n_m)
    !   all_trans(MORD, MORD, uhalf-1, n_m)
    real(c_double), pointer :: all_innov(:,:,:,:) => null()
    real(c_double), pointer :: all_trans(:,:,:,:) => null()
    integer :: win_all_innov = MPI_WIN_NULL
    integer :: win_all_trans = MPI_WIN_NULL

    ! cos(theta) at each upper-hemisphere ring (uhalf values, equator first)
    real(c_double), allocatable :: z_pts(:)

    ! MPI context stored at build time for reuse in draw_sphere_realisation.
    ! node_comm is MPI_COMM_NULL when running serially (no comm supplied).
    integer :: node_comm   = MPI_COMM_NULL
    integer :: node_rank   = 0
    integer :: node_npes   = 1
    ! m-mode strip owned by this rank: [j_lo, j_hi]  (1-based, inclusive)
    integer :: j_lo        = 1
    integer :: j_hi        = 0

    ! Ring strip for IRFFT (set at build time, reused in draw)
    integer :: row_lo      = 1
    integer :: row_hi      = 0
    integer :: n_row_local = 0

#ifdef MKL_AVAILABLE
    ! MKL DFTI descriptor for the batched 1-D real backward FFT.
    ! One descriptor per rank, created once in build_sphere_filter.
    ! Transform: folded complex CCE spectrum (n_fft) -> real (nphi),
    ! batched over n_row_local rings.
    type(DFTI_DESCRIPTOR), pointer :: dfti_irfft => null()
#endif

    logical :: initialized = .false.

  end type sphere_filter_type

  type, extends(abstract_random_fields) :: sphere_random_fields
    private
    type(sphere_filter_type) :: filter
    real(c_double), public :: c0 = 1.0d0
    real(c_double), public :: c2 = 0.0d0
    integer, public :: nz_sphere = 0
    integer, public :: nphi_sphere = 0
  contains
     procedure, public :: generate_2d_Random_field
     procedure, public :: generate_white_field
    procedure, public :: generate_2d_Random_field_seeded
    procedure, public :: finalize
    procedure, public :: initialized
  end type sphere_random_fields

  ! for example1_f in smerfs, so the random numbers are the same as those generated by python
  type :: mt19937_state
    integer(c_int32_t) :: mt(624)
    integer :: index = 625
  end type mt19937_state

  interface sphere_random_fields
    module procedure new_sphere_random_fields
  end interface sphere_random_fields

contains

  function sphere_random_fields_id(pert_param, pert_grid_f) result(id_string)
    type(pert_param_type), intent(in) :: pert_param
    type(grid_def_type), intent(in) :: pert_grid_f
    character(len=:), allocatable :: id_string
    character(len=256) :: id

    write(id,'(a,":",i0,":",i0,":",es16.8)') 'sphere', &
         pert_grid_f%N_lat, pert_grid_f%N_lon, pert_param%xcorr
    id_string = trim(id)
  end function sphere_random_fields_id

  function new_sphere_random_fields(pert_param, pert_grid_f, comm, rc) result(rf)
    type(sphere_random_fields) :: rf
    type(pert_param_type), intent(in) :: pert_param
    type(grid_def_type), intent(in) :: pert_grid_f
    integer, optional, intent(in) :: comm
    integer, optional, intent(out) :: rc
    real(c_double) :: xcorr_deg
    real(c_double) :: c2_use
    character(len=:), allocatable :: id_string
    integer :: comm_rank, mpierr

    if (present(rc)) rc = 0
    xcorr_deg = real(pert_param%xcorr, c_double)
    ! Sphere perturbations always use the fine grid; no coarsening is applied.
    rf%nphi_sphere = pert_grid_f%N_lon
    rf%nz_sphere = pert_grid_f%N_lat
    id_string = sphere_random_fields_id(pert_param, pert_grid_f)
    rf%ID_string = id_string
    if (xcorr_deg <= 0.0d0 .and. pert_param%ycorr <= 0.0) then
      ! land_pert calls generate_white_field for this case. No filter is needed.
      rf%filter%nphi = rf%nphi_sphere
      rf%filter%nz = rf%nz_sphere
      rf%filter%initialized = .true.
      return
    end if
    if (xcorr_deg <= 0.0d0) error stop 'sphere requires positive xcorr for correlated fields'
    if (rf%nphi_sphere < 4 .or. mod(rf%nphi_sphere,2) /= 0 .or. rf%nz_sphere < 2) &
      error stop 'sphere requires even N_lon >= 4 and N_lat >= 2'


    if (present(comm)) then
      call MPI_Comm_rank(comm, comm_rank, mpierr)
    else
      comm_rank = 0
    end if
    if (comm_rank == 0 .and. xcorr_deg > 0.0d0) then
      call xcorr_deg_to_c2(xcorr_deg, rf%nz_sphere, c2_use)
    else if (comm_rank == 0) then
      c2_use = 1.0d-12
    end if
    if (present(comm)) then
      call MPI_Bcast(c2_use, 1, MPI_DOUBLE_PRECISION, 0, comm, mpierr)
    end if
    rf%c2 = c2_use

    if (present(comm)) then
      call build_sphere_filter(rf%filter, rf%nz_sphere, rf%nphi_sphere, rf%c0, rf%c2, comm)
    else
      call build_sphere_filter(rf%filter, rf%nz_sphere, rf%nphi_sphere, rf%c0, rf%c2)
    end if
  end function new_sphere_random_fields

  subroutine generate_2d_Random_field(this, rseed, rfield, rfield2, lx, ly, dx, dy)
    class(sphere_random_fields), intent(inout) :: this
    integer, intent(inout) :: rseed(NRANDSEED)
    real, intent(out) :: rfield(:,:), rfield2(:,:)
    real, intent(in) :: lx, ly, dx, dy
    real(c_double), allocatable :: field_d(:,:), field2_d(:,:)
    real, allocatable :: output_field(:,:), output_field2(:,:)
    allocate(field_d(this%filter%nphi, this%filter%nz))
    allocate(field2_d(this%filter%nphi, this%filter%nz))
    call generate_2d_random_field_sphere(this%filter, rseed, field_d, field2_d, lx, ly, dx, dy)
    allocate(output_field(size(rfield,1), size(rfield,2)))
    allocate(output_field2(size(rfield2,1), size(rfield2,2)))
    call sphere_to_latlon_field(field_d, this%filter, size(rfield,1), size(rfield,2), output_field)
    call sphere_to_latlon_field(field2_d, this%filter, size(rfield2,1), size(rfield2,2), output_field2)
    rfield = output_field
    rfield2 = output_field2
    call normalize_unit_variance(rfield)
    call normalize_unit_variance(rfield2)
    deallocate(output_field, output_field2)
    deallocate(field_d, field2_d)
  end subroutine generate_2d_Random_field
  subroutine generate_white_field(this, rseed, rfield)
    class(sphere_random_fields), intent(inout) :: this
    integer, intent(inout) :: rseed(NRANDSEED)
    real, target, intent(out) :: rfield(:,:)
    real :: values(2)
    integer :: i, j

    do j = 1, size(rfield, 2)
      do i = 1, size(rfield, 1), 2
        call nr_gasdev(rseed, values)
        rfield(i,j) = values(1)
        if (i + 1 <= size(rfield, 1)) then
          rfield(i + 1,j) = values(2)
        end if
      end do
    end do
  end subroutine generate_white_field

  subroutine generate_2d_Random_field_seeded(this, seed, rfield, rfield2, lx, ly, dx, dy)
    class(sphere_random_fields), intent(inout) :: this
    integer, intent(in) :: seed
    real, intent(out) :: rfield(:,:), rfield2(:,:)
    real, intent(in) :: lx, ly, dx, dy
    ! The Python reference uses numpy.random.RandomState(seed), whose
    ! generator is MT19937.  Use the same raw 32-bit stream here instead of
    ! the project L'Ecuyer stream used by the general perturbation interface.
    call generate_2d_random_field_mt(this, seed, rfield, rfield2, lx, ly, dx, dy)
  end subroutine generate_2d_Random_field_seeded

  subroutine finalize(this, rc)
    class(sphere_random_fields), intent(inout) :: this
    integer, optional, intent(out) :: rc
    call finalize_sphere_filter(this%filter)
    if (present(rc)) rc = 0
  end subroutine finalize

  logical function initialized(this)
    class(sphere_random_fields), intent(in) :: this
    initialized = this%filter%initialized
  end function initialized

  ! ======================================================================
  !
  ! build_sphere_filter
  !
  ! Pre-compute state-space filter coefficients for the isotropic GRF
  ! with power spectrum  C_l = 1 / (c0 + c2*(l*(l+1))^2).
  !
  ! Parallel strategy:
  !   When comm is supplied it is split into node-local communicators, and
  !   node_comm is stored in sf for reuse in draw_sphere_realisation.  Each
  !   rank then runs the Riccati covariance / state-space
  !   chain only for its own m-mode strip and writes directly into the
  !   shared all_innov / all_trans windows. Integration is not duplicated
  !   within a node. A fence makes the complete filter
  !   visible to every rank.
  !
  ! Arguments:
  !   sf    (out)          : sphere_filter_type to be initialised
  !   nz    (in)           : number of iso-latitude rings (even)
  !   nphi  (in)           : number of phi points per ring (even)
  !   c0    (in)           : constant term of power spectrum denominator
  !   c2    (in)           : quartic term of power spectrum denominator
  !   comm  (in, optional) : MPI communicator (global, e.g. mpicomm)

  subroutine build_sphere_filter(sf, nz, nphi, c0, c2, comm)

    type(sphere_filter_type), intent(out) :: sf
    integer,        intent(in) :: nz, nphi
    real(c_double), intent(in) :: c0, c2
    integer, optional, intent(in) :: comm

    ! --- local variables ------------------------------------------------

    integer :: uhalf, n_m, m_cut, n_fft
    integer :: i, imode, rc
    integer :: mpierr
    integer :: n_m_build
    integer :: disp_unit
    integer(MPI_ADDRESS_KIND) :: win_size
    type(c_ptr) :: base_innov, base_trans
    real(c_double) :: dummy_real
    real(c_double) :: profile_start, profile_phase
    logical :: profile

    ! The new initializer returns the same 2-component storage layout.
    ! State component 2 is a scaled flux derivative, not an unscaled dg/dz.
    real(c_double) :: theta_val

    ! Temporary per-mode state-space arrays
    real(c_double), allocatable :: innov_m(:,:,:)   ! (MORD,MORD,uhalf)
    real(c_double), allocatable :: trans_m(:,:,:)   ! (MORD,MORD,uhalf-1)

    ! --------------------------------------------------------------------

    uhalf = nz/2 + 1
    n_fft = nphi/2 + 1
    call spherical_mode_cutoff(c0, c2, MODE_TAIL_TOL, m_cut, rc)
    if (rc /= 0) error stop 'invalid spherical mode cutoff'
    n_m = max(n_fft, m_cut+1)
    profile = sphere_profile_enabled()
    profile_start = MPI_Wtime()

    ! Store grid info in sf
    sf%nz    = nz
    sf%nphi  = nphi
    sf%n_m   = n_m
    sf%uhalf = uhalf
    sf%c0    = c0
    sf%c2    = c2

    ! --- MPI context ----------------------------------------------------
    ! Split the supplied communicator into node-local communicators.
    ! Store node_comm in sf so draw_sphere_realisation can reuse it.
    ! Compute the m-mode strip owned by this rank.

    if (present(comm)) then
      call MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, &
                               MPI_INFO_NULL, sf%node_comm, mpierr)
      call MPI_Comm_rank(sf%node_comm, sf%node_rank, mpierr)
      call MPI_Comm_size(sf%node_comm, sf%node_npes, mpierr)
    else
      sf%node_comm = MPI_COMM_NULL
      sf%node_rank = 0
      sf%node_npes = 1
    end if

    ! m-mode strip: distribute n_m modes as evenly as possible.  The same
    ! strip is used for the filter construction below and for the Kalman
    ! walk in draw_sphere_realisation.
    call compute_strip(n_m, sf%node_npes, sf%node_rank, sf%j_lo, sf%j_hi)
    n_m_build = max(0, sf%j_hi - sf%j_lo + 1)

    ! --- Allocate filter arrays (every rank needs the full tables) ------
    if (sf%node_comm /= MPI_COMM_NULL) then
      disp_unit = int(c_sizeof(dummy_real))

      win_size = 0_MPI_ADDRESS_KIND
      if (sf%node_rank == 0) win_size = int(MORD, MPI_ADDRESS_KIND) * &
           int(MORD, MPI_ADDRESS_KIND) * int(uhalf, MPI_ADDRESS_KIND) * &
           int(n_m, MPI_ADDRESS_KIND) * int(c_sizeof(dummy_real), MPI_ADDRESS_KIND)
      call MPI_Win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, &
                                   sf%node_comm, base_innov, sf%win_all_innov, mpierr)
      if (sf%node_rank /= 0) &
        call MPI_Win_shared_query(sf%win_all_innov, 0, win_size, disp_unit, base_innov, mpierr)
      call c_f_pointer(base_innov, sf%all_innov, [MORD, MORD, uhalf, n_m])

      win_size = 0_MPI_ADDRESS_KIND
      if (sf%node_rank == 0) win_size = int(MORD, MPI_ADDRESS_KIND) * &
           int(MORD, MPI_ADDRESS_KIND) * int(uhalf-1, MPI_ADDRESS_KIND) * &
           int(n_m, MPI_ADDRESS_KIND) * int(c_sizeof(dummy_real), MPI_ADDRESS_KIND)
      call MPI_Win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, &
                                   sf%node_comm, base_trans, sf%win_all_trans, mpierr)
      if (sf%node_rank /= 0) &
        call MPI_Win_shared_query(sf%win_all_trans, 0, win_size, disp_unit, base_trans, mpierr)
      call c_f_pointer(base_trans, sf%all_trans, [MORD, MORD, uhalf-1, n_m])
    else
      allocate(sf%all_innov(MORD, MORD, uhalf,   n_m))
      allocate(sf%all_trans(MORD, MORD, uhalf-1, n_m))
    end if

    ! --- z-points -------------------------------------------------------
    allocate(sf%z_pts(uhalf))
    do i = 1, uhalf
      theta_val   = (real(uhalf - i, c_double) + 0.5d0) * &
                    (PI / real(nz, c_double))
      sf%z_pts(i) = cos(theta_val)
    end do
    profile_phase = MPI_Wtime()
    if (profile) call sphere_profile_report('build setup and windows', &
                                            profile_phase-profile_start, comm)
  
    ! --- Parallel filter construction -----------------------------------
    ! Each rank builds the covariance/state-space chain only for its own
    ! m-mode strip and writes into the disjoint slices j_lo:j_hi of the
    ! shared filter tables.  The closing fence publishes the full tables.
    block
      if (sf%node_comm /= MPI_COMM_NULL) then
        if (sf%node_rank == 0) then
          sf%all_innov = 0.0d0
          sf%all_trans = 0.0d0
        end if
        call MPI_Win_fence(0, sf%win_all_innov, mpierr)
        call MPI_Win_fence(0, sf%win_all_trans, mpierr)
      else
        sf%all_innov = 0.0d0
        sf%all_trans = 0.0d0
      end if

      if (n_m_build > 0) then
        allocate(innov_m(MORD, MORD, uhalf))
        allocate(trans_m(MORD, MORD, uhalf-1))
        do imode = sf%j_lo - 1, sf%j_hi - 1
          call riccati_mode_filter(sf%z_pts, imode, c0, c2, innov_m, trans_m, rc)
          if (rc /= 0) then
            write(*,'(a,i0,a,i0)') &
              'build_sphere_filter: Riccati initializer failed at m=', imode, ', rc=', rc
            if (present(comm)) call MPI_Abort(comm, 1, mpierr)
            error stop 1
          end if
          sf%all_innov(:,:,:,imode+1) = innov_m
          sf%all_trans(:,:,:,imode+1) = trans_m
        end do
        deallocate(innov_m, trans_m)
      end if

      if (sf%node_comm /= MPI_COMM_NULL) then
        call MPI_Win_fence(0, sf%win_all_innov, mpierr)
        call MPI_Win_fence(0, sf%win_all_trans, mpierr)
      end if
    end block

    ! --- Pre-compute the ring strip owned by this rank for IRFFT --------
    ! Store in sf so draw_sphere_realisation can reuse it without
    ! recomputing compute_strip every call.
    call compute_strip(nz, sf%node_npes, sf%node_rank, sf%row_lo, sf%row_hi)
    sf%n_row_local = sf%row_hi - sf%row_lo + 1

#ifdef MKL_AVAILABLE
    ! --- Create MKL DFTI descriptor for batched 1-D real IRFFT ----------
    ! Each rank transforms sf%n_row_local rings of length nphi (single prec).
    !
    ! Transform: DFTI_REAL domain, DFTI_BACKWARD (complex half-spec -> real)
    ! Storage:   in-place (the default), so input and output share one
    !            buffer.  Each ring occupies 2*n_m real(c_float) elements:
    !            MKL reads the n_m complex CCE values and overwrites them
    !            with the nphi real output values in place.
    !
    ! Batching over n_row_local rings:
    !   DFTI_NUMBER_OF_TRANSFORMS = n_row_local
    !   DFTI_INPUT_DISTANCE  = n_m    (complex elements between ring starts)
    !   DFTI_OUTPUT_DISTANCE = 2*n_m  (real    elements between ring starts)
    !
    ! DFTI_COMPLEX_STORAGE defaults to DFTI_COMPLEX_COMPLEX for real domain
    ! backward transforms, meaning the input is treated as n_m complex values
    ! (the CCE half-spectrum).  No explicit storage-format call is needed.
    block
      integer :: dfti_stat
      if (sf%n_row_local > 0) then
        dfti_stat = DftiCreateDescriptor(sf%dfti_irfft, &
                                         DFTI_SINGLE, DFTI_REAL, 1, nphi)
        if (dfti_stat /= DFTI_NO_ERROR) then
          write(*,'(a,i0)') &
            'build_sphere_filter: DftiCreateDescriptor failed, stat=', dfti_stat
          stop 1
        end if
        dfti_stat = DftiSetValue(sf%dfti_irfft, &
                                 DFTI_NUMBER_OF_TRANSFORMS, sf%n_row_local)
        dfti_stat = DftiSetValue(sf%dfti_irfft, DFTI_INPUT_DISTANCE,  n_fft)
        dfti_stat = DftiSetValue(sf%dfti_irfft, DFTI_OUTPUT_DISTANCE, 2*n_fft)
        dfti_stat = DftiCommitDescriptor(sf%dfti_irfft)
        if (dfti_stat /= DFTI_NO_ERROR) then
          write(*,'(a,i0)') &
            'build_sphere_filter: DftiCommitDescriptor failed, stat=', dfti_stat
          stop 1
        end if
      end if
    end block
#endif

    sf%initialized = .true.
    if (profile) call sphere_profile_report('build total', MPI_Wtime()-profile_start, comm)

  end subroutine build_sphere_filter

  ! ======================================================================
  !
  ! draw_sphere_realisation
  !
  ! Generate one realisation of the isotropic GRF on the sphere.
  !
  ! Parallel strategy (m-mode distribution):
  !
  !   The n_m Fourier m-modes are split evenly across the node_npes ranks
  !   of the node communicator stored in sf (established at build time).
  !   Each rank:
  !     1. Generates noise only for its m-mode strip j_lo:j_hi.
  !     2. Runs the Kalman walk independently for those m-modes.
  !     3. Writes res_cplx(j_lo:j_hi, :) into a shared MPI window.
  !   After MPI_Win_fence, all ranks see the complete res_cplx(n_m, nz).
  !   The IRFFT rows (nz rings) are then distributed across ranks:
  !     Each rank transforms rfield_out(:, row_lo:row_hi) into a shared
  !     rfield window.  After a second fence every rank holds the full
  !     rfield_out(nphi, nz).
  !
  !   When sf%node_comm == MPI_COMM_NULL (serial build / no comm), the
  !   subroutine runs the same walk over a single strip covering all
  !   modes, using private allocatables instead of shared windows.
  !
  ! Arguments:
  !   sf           (in)  : initialised sphere_filter_type
  !   rand_ints    (in)  : array of random 32-bit integers.  Every rank
  !                        draws from the same global stream and selects
  !                        its own strip, so the size must be >= the
  !                        global requirement nz*MORD*n_m*2 irrespective
  !                        of the number of ranks.
  !   n_rand_got   (in)  : length of rand_ints supplied
  !   rfield_out   (out) : real(c_double) field (sf%nphi, sf%nz)
  !                        On return every rank holds the complete field.

  subroutine draw_sphere_realisation(sf, rand_ints, n_rand_got, rfield_out, rseed)

    type(sphere_filter_type),            intent(in)  :: sf
    integer,                             intent(in)  :: n_rand_got
    integer(c_int32_t),                  intent(in)  :: rand_ints(n_rand_got)
    real(c_double),                      intent(out) :: rfield_out(sf%nphi, sf%nz)
    integer, optional, intent(inout) :: rseed(NRANDSEED)
    real :: gaussian_pair(2)

    ! --- local variables ------------------------------------------------

    integer :: uhalf, n_m, nz, nphi, n_fft
    integer :: j_lo, j_hi, n_m_local           ! m-mode strip for this rank
    integer :: row_lo, row_hi, n_row_local     ! ring strip for IRFFT (from sf)
    integer :: n_rand_needed_global
    integer :: i, j, ip, iq, nn, rc
    integer :: skip
    logical :: parallel

    real(c_float),    allocatable :: noise_real(:)       ! standard normals
    complex(c_float), allocatable :: noise(:,:,:)        ! (MORD, n_m_local, nz)
    complex(c_float), allocatable :: fp_f(:,:)           ! (MORD, n_m_local)
    complex(c_float), allocatable :: fp_fstart(:,:)      ! (MORD, n_m_local)
    complex(c_float), allocatable :: fp_fprev(:,:)       ! evolving state

    ! Half-spectrum and real-space field (shared windows in parallel,
    ! private allocatables in serial)
    complex(c_float), allocatable, target :: res_cplx_priv(:,:)    ! (n_m, nz)  serial
    real(c_double),   allocatable, target :: rfield_priv(:,:)       ! (nphi, nz) serial

    complex(c_float), pointer :: res_cplx(:,:)   ! points into window or priv
    real(c_double),   pointer :: rfield_win(:,:)  ! points into window or priv

    ! MPI shared-window variables
    integer :: win_res, win_rf
    integer :: mpierr
    integer(MPI_ADDRESS_KIND) :: win_size
    integer :: disp_unit
    type(c_ptr) :: base_res, base_rf
    complex(c_float) :: dummy_c
    real(c_double)   :: dummy_r
    real(c_double) :: profile_start, profile_phase
    logical :: profile

    ! --------------------------------------------------------------------

    uhalf = sf%uhalf
    n_m   = sf%n_m
    nz    = sf%nz
    nphi  = sf%nphi
    n_fft = nphi/2 + 1
    j_lo  = sf%j_lo
    j_hi  = sf%j_hi

    parallel = (sf%node_comm /= MPI_COMM_NULL)
    profile = sphere_profile_enabled()
    profile_start = MPI_Wtime()

    n_m_local = j_hi - j_lo + 1
    if (n_m_local <= 0) n_m_local = 0   ! edge case: more ranks than modes

    ! Ring strip pre-computed at build time; reuse from sf
    row_lo      = sf%row_lo
    row_hi      = sf%row_hi
    n_row_local = sf%n_row_local

    ! --- 1. Allocate or map shared windows for res_cplx and rfield_out --

    if (parallel) then
      ! --- res_cplx shared window: shape (n_m, nz), complex(c_float) ---
      disp_unit = 4
      win_size  = 0_MPI_ADDRESS_KIND
      if (sf%node_rank == 0) &
        win_size = int(n_m, MPI_ADDRESS_KIND) * int(nz, MPI_ADDRESS_KIND) &
                   * int(c_sizeof(dummy_c), MPI_ADDRESS_KIND)
      call MPI_Win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, &
                                   sf%node_comm, base_res, win_res, mpierr)
      if (sf%node_rank /= 0) &
        call MPI_Win_shared_query(win_res, 0, win_size, disp_unit, base_res, mpierr)
      call MPI_Win_fence(0, win_res, mpierr)
      call c_f_pointer(base_res, res_cplx, [n_m, nz])

      ! --- rfield_out shared window: shape (nphi, nz), real(c_double) ---
      win_size = 0_MPI_ADDRESS_KIND
      if (sf%node_rank == 0) &
        win_size = int(nphi, MPI_ADDRESS_KIND) * int(nz, MPI_ADDRESS_KIND) &
                   * int(c_sizeof(dummy_r), MPI_ADDRESS_KIND)
      call MPI_Win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, &
                                   sf%node_comm, base_rf, win_rf, mpierr)
      if (sf%node_rank /= 0) &
        call MPI_Win_shared_query(win_rf, 0, win_size, disp_unit, base_rf, mpierr)
      call MPI_Win_fence(0, win_rf, mpierr)
      call c_f_pointer(base_rf, rfield_win, [nphi, nz])
    else
      ! Serial: use private allocatables and point to them.  The strips
      ! stored in sf already cover everything (node_npes == 1).
      allocate(res_cplx_priv(n_m, nz))
      allocate(rfield_priv  (nphi, nz))
      res_cplx => res_cplx_priv
      rfield_win => rfield_priv
    end if
    profile_phase = MPI_Wtime()
    if (profile) call sphere_profile_report('draw windows', profile_phase-profile_start, sf%node_comm)

    ! --- 2. Noise generation for the local m-mode strip -----------------
    ! The whole global stream is drawn on every rank so that each rank can
    ! index its own strip out of it and the result is independent of the
    ! rank count.  The south-hemisphere walk reuses noise at indices
    ! uhalf+1..nz, so the full nz rings of noise are required.
    n_rand_needed_global = nz * MORD * n_m * 2

    allocate(noise_real(n_rand_needed_global))
    allocate(noise    (MORD, n_m_local, nz))
    allocate(fp_f     (MORD, n_m_local))
    allocate(fp_fstart(MORD, n_m_local))
    allocate(fp_fprev (MORD, n_m_local))

    profile_phase = MPI_Wtime()
    ! Production: use the existing restartable NR Gaussian generator directly.
    ! Draw the same global stream on every rank, including ranks with empty strips.
    if (present(rseed)) then
      do nn = 1, n_rand_needed_global, 2
        call nr_gasdev(rseed, gaussian_pair)
        noise_real(nn:nn+1) = real(gaussian_pair, c_float)
      end do
    else
      ! Raw full-width MT words are retained only for the benchmark interface.
      call zigg_f(n_rand_needed_global, n_rand_got, rand_ints, noise_real, rc)
      if (rc /= 0) then
        write(*,'(a,i0)') 'draw_sphere_realisation: Gaussian stream exhausted, rc=', rc
        if (parallel) call MPI_Abort(sf%node_comm, 1, mpierr)
        error stop 1
      end if
    end if

    if (n_m_local > 0) then
      ! Pack into complex noise: noise(ip, jlocal, i) = (re + i*im) / sqrt(2)
      ! Loop over all nz rings so south-hemisphere indices (uhalf+1..nz)
      ! are filled from the same noise stream as the reference Python code.
      do i = 1, nz
        do ip = 1, MORD
          do j = 1, n_m_local
            nn = (i - 1) * MORD * n_m + (ip - 1) * n_m + (j_lo + j - 2) + 1
            noise(ip, j, i) = cmplx(noise_real(2*nn-1), noise_real(2*nn), c_float) &
                               * real(0.5d0**0.5d0, c_float)
          end do
        end do
      end do

      ! --- 3. Kalman walk for local m-modes j_lo:j_hi -------------------
      ! Upward walk: equator (ring 1 in sf%z_pts) -> pole (ring uhalf)

      ! Initialise at ring 1
      do j = 1, n_m_local
        do ip = 1, MORD
          fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
          do iq = 1, MORD
            fp_f(ip, j) = fp_f(ip, j) &
               + real(sf%all_innov(iq, ip, 1, j_lo+j-1), c_float) * noise(iq, j, 1)
          end do
        end do
      end do
      fp_fstart = fp_f
      fp_fprev = fp_f
      res_cplx(j_lo:j_hi, uhalf) = fp_f(1, :)

      do i = 1, uhalf-1
        do j = 1, n_m_local
          do ip = 1, MORD
            fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
            do iq = 1, MORD
              fp_f(ip, j) = fp_f(ip, j) &
                 + real(sf%all_trans(iq, ip, i,   j_lo+j-1), c_float) * fp_fprev(iq, j) &
                 + real(sf%all_innov(iq, ip, i+1, j_lo+j-1), c_float) * noise(iq, j, i+1)
            end do
          end do
        end do
        fp_fprev = fp_f
        res_cplx(j_lo:j_hi, uhalf - i) = fp_f(1, :)
      end do

      ! Reflect across the equator
      do ip = 1, MORD
        if (mod(ip-1, 2) == 0) then
          fp_f(ip, :) = fp_fstart(ip, :)
        else
          fp_f(ip, :) = -fp_fstart(ip, :)
        end if
      end do
      fp_fprev = fp_f
      skip = 1 - mod(nz, 2)

      do i = 0, nz - uhalf - 1
        do j = 1, n_m_local
          do ip = 1, MORD
            fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
            do iq = 1, MORD
              fp_f(ip, j) = fp_f(ip, j) &
                  + real(sf%all_trans(iq, ip, i+skip+1, j_lo+j-1), c_float) * fp_fprev(iq, j) &
                  + real(sf%all_innov(iq, ip, i+skip+2, j_lo+j-1), c_float) * noise(iq, j, uhalf+i+1)
            end do
          end do
        end do
        fp_fprev = fp_f
        res_cplx(j_lo:j_hi, uhalf + i + 1) = fp_f(1, :)
      end do

      ! Scale DC mode (m=0, j_lo=1 only on the rank that owns mode 0)
      if (j_lo == 1) then
        do i = 1, nz
          res_cplx(1, i) = cmplx(real(res_cplx(1, i)) * sqrt(2.0), 0.0, c_float)
        end do
      end if

    end if  ! n_m_local > 0

    if (profile) call sphere_profile_report('draw noise, packing, and spectral walk', &
                                            MPI_Wtime()-profile_phase, sf%node_comm)
    profile_phase = MPI_Wtime()

    ! --- 4. Fence: all ranks have written their m-mode strip ------------
    if (parallel) then
      call MPI_Win_fence(0, win_res, mpierr)
    end if
    if (profile) call sphere_profile_report('draw spectrum fence', &
                                            MPI_Wtime()-profile_phase, sf%node_comm)
    profile_phase = MPI_Wtime()

    ! --- 5. Distributed IRFFT: each rank transforms its ring strip ------
    ! MKL path: batched 1-D real backward FFT using the pre-built descriptor.
    ! Fallback: pure-Fortran O(nphi^2) DFT (NR path or n_row_local==0).
#ifdef MKL_AVAILABLE
    if (n_row_local > 0) then
      block
        ! In-place real backward FFT using the pre-built descriptor.
        ! The work buffer holds n_row_local rings, each padded to 2*n_m
        ! real(c_float) elements (= n_m complex elements).
        ! MKL reads the folded half-spectrum from the first n_fft complex
        ! positions (= 2*n_fft reals) and writes the real output (nphi reals)
        ! starting at the same address.
        real(c_float), allocatable, target :: work(:)   ! (2*n_fft * n_row_local)
        type(c_ptr) :: cptr_w
        complex(c_float), pointer :: cplx_view(:,:)  ! (n_fft, n_row_local)
        real(c_float),    pointer :: real_view(:,:)   ! (2*n_fft, n_row_local)
        integer :: dfti_stat, row

        allocate(work(2 * n_fft * n_row_local))

        ! Map complex view onto work buffer, copy half-spectrum in
        cptr_w = c_loc(work(1))
        call c_f_pointer(cptr_w, cplx_view, [n_fft, n_row_local])
        call fold_fourier_modes(res_cplx(:, row_lo:row_hi), nphi, cplx_view)

        ! In-place backward FFT
        dfti_stat = DftiComputeBackward(sf%dfti_irfft, work)
        if (dfti_stat /= DFTI_NO_ERROR) then
          write(*,'(a,i0)') &
            'draw_sphere_realisation: DftiComputeBackward failed, stat=', dfti_stat
          stop 1
        end if

        ! Real output occupies first nphi reals of each padded row (stride 2*n_fft)
        call c_f_pointer(cptr_w, real_view, [2*n_fft, n_row_local])
        do row = 1, n_row_local
          rfield_win(:, row_lo + row - 1) = real(real_view(1:nphi, row), c_double)
        end do

        deallocate(work)
      end block
    end if
#else
    if (n_row_local > 0) then
      block
        complex(c_float), allocatable :: folded(:,:)
        allocate(folded(n_fft,n_row_local))
        call fold_fourier_modes(res_cplx(:,row_lo:row_hi), nphi, folded)
        call irfft_rows(folded, n_row_local, n_fft, nphi, rfield_win(:,row_lo:row_hi))
        deallocate(folded)
      end block
    end if
#endif

    ! --- 6. Fence: all ranks have written their ring strip ---------------
    if (parallel) call MPI_Win_fence(0, win_rf, mpierr)
    if (profile) call sphere_profile_report('draw inverse FFT and fence', &
                                            MPI_Wtime()-profile_phase, sf%node_comm)
    profile_phase = MPI_Wtime()

    ! --- 7. Normalise to unit variance (global empirical) ----------------
    ! The IRFFT output has a variance that depends on nz, nphi, MORD, and
    ! the power spectrum.  We normalise the full field to zero mean and unit
    ! variance empirically so that the std passed from land_pert.F90 can be
    ! applied directly.
    !
    ! IMPORTANT: this normalisation is computed and applied on the FULL
    ! rfield_win (all nz rings) AFTER the fence, so all rings are scaled by
    ! the same factor and the relative ring-to-ring structure is preserved.
    ! Each rank computes the sum-of-squares over its own IRFFT ring strip
    ! (row_lo:row_hi); MPI_Allreduce gathers the global total.  Because
    ! every rank holds the full rfield_win after the fence, each rank then
    ! divides its strip by the global norm_fac.  A second fence ensures that
    ! rfield_win is fully normalised before being copied to rfield_out.
    block
      real(c_double) :: local_sum, global_sum, global_mean
      real(c_double) :: local_sum2, global_sum2, global_var, norm_fac
      integer :: ntot, mpierr2

      if (n_row_local > 0) then
        local_sum = sum(rfield_win(:, row_lo:row_hi))
      else
        local_sum = 0.0d0
      end if
      ntot = nphi * nz

      if (parallel) then
        call MPI_Allreduce(local_sum, global_sum, 1, &
                           MPI_DOUBLE_PRECISION, MPI_SUM, sf%node_comm, mpierr2)
      else
        global_sum = local_sum
      end if
      global_mean = global_sum / real(ntot, c_double)

      if (n_row_local > 0) &
        rfield_win(:, row_lo:row_hi) = rfield_win(:, row_lo:row_hi) - global_mean

      if (n_row_local > 0) then
        local_sum2 = sum(rfield_win(:, row_lo:row_hi)**2)
      else
        local_sum2 = 0.0d0
      end if

      if (parallel) then
        call MPI_Allreduce(local_sum2, global_sum2, 1, &
                           MPI_DOUBLE_PRECISION, MPI_SUM, sf%node_comm, mpierr2)
      else
        global_sum2 = local_sum2
      end if

      global_var = global_sum2 / real(ntot, c_double)
      norm_fac   = sqrt(global_var)

      if (n_row_local > 0 .and. norm_fac > 0.0d0) &
        rfield_win(:, row_lo:row_hi) = rfield_win(:, row_lo:row_hi) / norm_fac
    end block
    if (profile) call sphere_profile_report('draw normalization', &
                                            MPI_Wtime()-profile_phase, sf%node_comm)

    ! --- 8. Final fence + copy to caller --------------------------------
    if (parallel) then
      call MPI_Win_fence(0, win_rf, mpierr)
      rfield_out = rfield_win
      call MPI_Win_free(win_res, mpierr)
      call MPI_Win_free(win_rf,  mpierr)
    else
      rfield_out = rfield_win
      deallocate(res_cplx_priv, rfield_priv)
      nullify(res_cplx, rfield_win)
    end if

    deallocate(noise_real, noise, fp_f, fp_fstart, fp_fprev)
    if (profile) call sphere_profile_report('draw total', MPI_Wtime()-profile_start, sf%node_comm)

  end subroutine draw_sphere_realisation

  ! ======================================================================
  !
  ! finalize_sphere_filter

  subroutine finalize_sphere_filter(sf)
    type(sphere_filter_type), intent(inout) :: sf
    integer :: mpierr

    if (sf%win_all_innov /= MPI_WIN_NULL) then
      call MPI_Win_free(sf%win_all_innov, mpierr)
      sf%win_all_innov = MPI_WIN_NULL
      nullify(sf%all_innov)
    else if (associated(sf%all_innov)) then
      deallocate(sf%all_innov)
    end if
    if (sf%win_all_trans /= MPI_WIN_NULL) then
      call MPI_Win_free(sf%win_all_trans, mpierr)
      sf%win_all_trans = MPI_WIN_NULL
      nullify(sf%all_trans)
    else if (associated(sf%all_trans)) then
      deallocate(sf%all_trans)
    end if
    if (allocated(sf%z_pts))     deallocate(sf%z_pts)
#ifdef MKL_AVAILABLE
    if (associated(sf%dfti_irfft)) then
      mpierr = DftiFreeDescriptor(sf%dfti_irfft)
      nullify(sf%dfti_irfft)
    end if
#endif
    if (sf%node_comm /= MPI_COMM_NULL) then
      call MPI_Comm_free(sf%node_comm, mpierr)
      sf%node_comm = MPI_COMM_NULL
    end if
    sf%initialized = .false.
  end subroutine finalize_sphere_filter

  ! ======================================================================
  !
  ! sphere_to_latlon_field
  !
  ! Remap sphere-native field rfield_sphere(nphi, nz) onto the standard
  ! equidistant lat-lon grid rfield_latlon(N_lon, N_lat).
  ! Bilinear interpolation in z=sin(latitude) and periodic longitude.

  subroutine sphere_to_latlon_field(rfield_sphere, sf, N_lon, N_lat, rfield_latlon)

    type(sphere_filter_type), intent(in)  :: sf
    real(c_double),           intent(in)  :: rfield_sphere(sf%nphi, sf%nz)
    integer,                  intent(in)  :: N_lon, N_lat
    real,                     intent(out) :: rfield_latlon(N_lon, N_lat)

    integer :: i_lon, j_lat, k1, k2, p1, p2, k
    real(c_double) :: lat_deg, lat_rad, z_target, z1, z2, wz
    real(c_double) :: phi_pos, wp
    real(c_double) :: v1, v2
    do j_lat = 1, N_lat
      lat_deg  = -90.0d0 + (real(j_lat, c_double) - 0.5d0) * 180.0d0 / real(N_lat, c_double)
      lat_rad  = lat_deg * MAPL_DEGREES_TO_RADIANS_R8
      z_target = sin(lat_rad)

      if (z_target >= ring_z(1, sf)) then
        k1 = 1; k2 = 1; wz = 0.0d0
      else if (z_target <= ring_z(sf%nz, sf)) then
        k1 = sf%nz; k2 = sf%nz; wz = 0.0d0
      else
        k1 = 1
        do k = 1, sf%nz - 1
          if (ring_z(k, sf) >= z_target .and. z_target >= ring_z(k+1, sf)) then
            k1 = k
            exit
          end if
        end do
        k2 = k1 + 1
        z1 = ring_z(k1, sf)
        z2 = ring_z(k2, sf)
        wz = (z1 - z_target) / (z1 - z2)
      end if

      do i_lon = 1, N_lon
        phi_pos = 0.5d0 + (real(i_lon, c_double) - 0.5d0) * real(sf%nphi, c_double) / real(N_lon, c_double)
        p1 = floor(phi_pos)
        wp = phi_pos - real(p1, c_double)
        p1 = modulo(p1 - 1, sf%nphi) + 1
        p2 = modulo(p1, sf%nphi) + 1
        v1 = (1.0d0 - wz) * rfield_sphere(p1,k1) + wz * rfield_sphere(p1,k2)
        v2 = (1.0d0 - wz) * rfield_sphere(p2,k1) + wz * rfield_sphere(p2,k2)
        rfield_latlon(i_lon, j_lat) = real((1.0d0-wp)*v1 + wp*v2)
      end do
    end do

  end subroutine sphere_to_latlon_field

  ! ======================================================================
  !
  ! xcorr_deg_to_c2
  !
  ! Derive SMERFS power-spectrum coefficient c2 from a user-supplied
  ! horizontal correlation length xcorr_deg [degrees].
  !
  ! Match the flat-FFT path's correlation LEVEL at the requested distance:
  ! C(xcorr_deg)/C(0) = exp(-0.5). The rational spherical spectrum is not
  ! a Gaussian distance kernel, so its full curve and e-folding distance
  ! must NOT be inferred from a Gaussian formula.
  ! Calibration uses a grid-independent Legendre sum with a tail bound.

  subroutine xcorr_deg_to_c2(xcorr_deg, nz, c2_out)
    real(c_double), intent(in) :: xcorr_deg
    integer, intent(in) :: nz  ! Retained for private-call compatibility, not used in calibration.
    real(c_double), intent(out) :: c2_out
    integer :: rc_local, mpierr
    call calibrate_spherical_c2(xcorr_deg, 1.0d0, c2_out, rc_local)
    if (rc_local /= 0) then
      write(*,'(a,f12.6,a,i0)') 'xcorr calibration failed: ', xcorr_deg, ', rc=', rc_local
      call MPI_Abort(MPI_COMM_WORLD, 1, mpierr)
      error stop 1
    end if
  end subroutine xcorr_deg_to_c2

  subroutine generate_2d_random_field_sphere(this, rseed, rfield, rfield2, lx, ly, dx, dy)
    type(sphere_filter_type), intent(in) :: this
    integer, intent(inout) :: rseed(NRANDSEED)
    real(c_double), intent(out) :: rfield(this%nphi, this%nz), rfield2(this%nphi, this%nz)
    real, intent(in) :: lx, ly, dx, dy
    integer(c_int32_t) :: unused_ints(1) = 0
    call draw_sphere_realisation(this, unused_ints, 1, rfield, rseed=rseed)
    call draw_sphere_realisation(this, unused_ints, 1, rfield2, rseed=rseed)
  end subroutine generate_2d_random_field_sphere

  ! Generate a field using the raw MT19937 stream used by numpy's
  ! RandomState(seed).  The returned integers have the same bit pattern as
  ! numpy's uint32 values and are consumed by the same C Ziggurat routine.
  subroutine generate_2d_random_field_mt(this, seed, rfield, rfield2, lx, ly, dx, dy)
    class(sphere_random_fields), intent(inout) :: this
    integer, intent(in) :: seed
    real, intent(out) :: rfield(:,:)
    real, intent(out) :: rfield2(:,:)
    real, intent(in) :: lx, ly, dx, dy
    integer(c_int32_t), allocatable :: ints(:)
    integer :: n_needed, n_got
    type(mt19937_state) :: state
    real(c_double), allocatable :: field_d(:,:), field2_d(:,:)

    n_needed = this%filter%nz * MORD * this%filter%n_m * 2
    n_got = n_needed + n_needed / 32 + 64
    allocate(ints(n_got))
    allocate(field_d(this%filter%nphi, this%filter%nz))
    allocate(field2_d(this%filter%nphi, this%filter%nz))
    call mt19937_init(state, seed)
    call mt19937_fill(state, ints)
    call draw_sphere_realisation(this%filter, ints, n_got, field_d)
    call mt19937_fill(state, ints)
    call draw_sphere_realisation(this%filter, ints, n_got, field2_d)
    rfield = real(field_d)
    rfield2 = real(field2_d)
    call normalize_unit_variance(rfield)
    call normalize_unit_variance(rfield2)
    deallocate(ints, field_d, field2_d)
  end subroutine generate_2d_random_field_mt

  subroutine mt19937_init(state, seed)
    type(mt19937_state), intent(out) :: state
    integer, intent(in) :: seed
    integer :: i
    integer(c_int64_t) :: value
    integer(c_int64_t), parameter :: modulus = 4294967296_c_int64_t

    state%mt(1) = int(seed, c_int32_t)
    do i = 2, 624
      value = 1812433253_c_int64_t * &
              iand(int(ieor(state%mt(i-1), ishft(state%mt(i-1), -30)), c_int64_t), &
                   modulus - 1_c_int64_t) + i - 1
      state%mt(i) = int(modulo(value, modulus), c_int32_t)
    end do
    state%index = 625
  end subroutine mt19937_init

  subroutine mt19937_twist(state)
    type(mt19937_state), intent(inout) :: state
    integer :: i, j
    integer(c_int32_t) :: y
    integer(c_int32_t), parameter :: upper_mask = int(z'80000000', c_int32_t)
    integer(c_int32_t), parameter :: lower_mask = int(z'7fffffff', c_int32_t)
    integer(c_int32_t), parameter :: matrix_a = int(z'9908b0df', c_int32_t)

    do i = 1, 624
      j = i + 1
      if (j > 624) j = 1
      y = ior(iand(state%mt(i), upper_mask), iand(state%mt(j), lower_mask))
      j = i + 397
      if (j > 624) j = j - 624
      state%mt(i) = ieor(state%mt(j), ishft(y, -1))
      if (iand(y, 1_c_int32_t) /= 0_c_int32_t) &
        state%mt(i) = ieor(state%mt(i), matrix_a)
    end do
    state%index = 1
  end subroutine mt19937_twist

  function mt19937_next(state) result(value)
    type(mt19937_state), intent(inout) :: state
    integer(c_int32_t) :: value
    integer(c_int32_t) :: y

    if (state%index > 624) call mt19937_twist(state)
    y = state%mt(state%index)
    state%index = state%index + 1
    y = ieor(y, iand(ishft(y, -11), int(z'ffffffff', c_int32_t)))
    y = ieor(y, iand(ishft(y, 7),  int(z'9d2c5680', c_int32_t)))
    y = ieor(y, iand(ishft(y, 15), int(z'efc60000', c_int32_t)))
    y = ieor(y, ishft(y, -18))
    value = y
  end function mt19937_next

  subroutine mt19937_fill(state, values)
    type(mt19937_state), intent(inout) :: state
    integer(c_int32_t), intent(out) :: values(:)
    integer :: i
    do i = 1, size(values)
      values(i) = mt19937_next(state)
    end do
  end subroutine mt19937_fill

  ! ======================================================================
  ! Private: normalize_unit_variance
  !
  ! Shift field to zero mean and scale to unit variance, so that the std
  ! passed down from land_pert.F90 can be applied directly.  A zero-norm
  ! field is left untouched.

  subroutine normalize_unit_variance(field)
    real, intent(inout) :: field(:,:)
    real :: field_mean, field_norm

    field_mean = sum(field) / real(size(field))
    field = field - field_mean
    field_norm = sqrt(sum(field**2) / real(size(field)))
    if (field_norm > 0.0) field = field / field_norm
  end subroutine normalize_unit_variance

  ! ======================================================================
  ! Private: compute_strip
  !
  ! Distribute n total items across npes ranks; return the 1-based
  ! inclusive range [lo, hi] owned by rank r (0-based).
  ! Remainder items go to the first (mod) ranks.

  subroutine compute_strip(n, npes, r, lo, hi)
    integer, intent(in)  :: n, npes, r
    integer, intent(out) :: lo, hi
    integer :: base, rem
    base = n / npes
    rem  = mod(n, npes)
    ! ranks 0..rem-1 get (base+1) items; ranks rem..npes-1 get base items
    if (r < rem) then
      lo = r * (base + 1) + 1
      hi = lo + base
    else
      lo = rem * (base + 1) + (r - rem) * base + 1
      hi = lo + base - 1
    end if
  end subroutine compute_strip

  logical function sphere_profile_enabled()
    character(len=8) :: value
    integer :: length, status

    call get_environment_variable('SPHERE_RF_PROFILE', value, length, status)
    sphere_profile_enabled = status == 0 .and. length > 0 .and. value(1:1) /= '0'
  end function sphere_profile_enabled

  subroutine sphere_profile_report(label, seconds, comm)
    character(*), intent(in) :: label
    real(c_double), intent(in) :: seconds
    integer, optional, intent(in) :: comm
    integer :: rank, mpierr
    real(c_double) :: max_seconds

    if (present(comm)) then
      call MPI_Reduce(seconds, max_seconds, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, mpierr)
      call MPI_Comm_rank(comm, rank, mpierr)
      if (rank /= 0) return
    else
      max_seconds = seconds
    end if
    write(*,'(a,a,a,f10.4,a)') 'SPHERE_PROFILE: ', trim(label), ' max=', max_seconds, ' s'
  end subroutine sphere_profile_report

  ! ======================================================================
  ! Private: native ring z coordinate

  real(c_double) function ring_z(k, sf)
    integer, intent(in) :: k
    type(sphere_filter_type), intent(in) :: sf

    if (k <= sf%uhalf) then
      ring_z = sf%z_pts(sf%uhalf + 1 - k)
    else
      ring_z = -sf%z_pts(k - sf%uhalf + 1 + (1 - mod(sf%nz,2)))
    end if
  end function ring_z

  ! Fold positive Fourier modes into the real-FFT bins. Each rank owns
  ! disjoint output rows, so alias additions require no shared-bin reduction.
  ! Implements g0 + 2*Re(sum_(m>=1) g_m*exp(i*m*phi)) at output longitudes.
  subroutine fold_fourier_modes(full, nphi, folded)
    integer, intent(in) :: nphi
    complex(c_float), intent(in) :: full(:,:)
    complex(c_float), intent(out) :: folded(nphi/2+1,size(full,2))
    integer :: row,m,k
    folded=cmplx(0.0,0.0,c_float)
    do row=1,size(full,2)
      folded(1,row)=cmplx(real(full(1,row)),0.0,c_float)
      do m=1,size(full,1)-1
        k=modulo(m,nphi)
        if(k==0) then
          folded(1,row)=folded(1,row)+cmplx(2.0*real(full(m+1,row)),0.0,c_float)
        else if(k==nphi/2) then
          folded(nphi/2+1,row)=folded(nphi/2+1,row)+cmplx(2.0*real(full(m+1,row)),0.0,c_float)
        else if(k<nphi/2) then
          folded(k+1,row)=folded(k+1,row)+full(m+1,row)
        else
          folded(nphi-k+1,row)=folded(nphi-k+1,row)+conjg(full(m+1,row))
        end if
      end do
    end do
  end subroutine fold_fourier_modes

  ! ======================================================================
  ! Private: irfft_rows
  !
  ! Inverse real FFT of each column (ring) of the half-spectrum.
  !   inp(nm, nrows) complex -> out(nph, nrows) real
  ! Pure Fortran O(nph^2) DFT.

  subroutine irfft_rows(inp, nrows, nm, nph, out)

    integer,          intent(in)  :: nrows, nm, nph
    complex(c_float), intent(in)  :: inp(nm, nrows)
    real(c_double),   intent(out) :: out(nph, nrows)

    integer        :: row, j, k
    real(c_double) :: angle, re_sum, two_pi_over_n
    complex(c_double_complex) :: Xk
    ! The SMERFS Kalman filter is calibrated to produce a unit-variance
    ! field when using the DFT convention without the 1/N normalisation
    ! factor (matching example1_f.F90 and the Python `res *= p` step).
    two_pi_over_n = 2.0d0 * PI / real(nph, c_double)

    do row = 1, nrows
      do j = 0, nph-1
        re_sum = real(inp(1, row), c_double)
        do k = 1, nm-2
          angle  = two_pi_over_n * real(j*k, c_double)
          Xk = cmplx(real(inp(k+1, row), c_double), &
                     real(aimag(inp(k+1, row)), c_double), c_double_complex)
          re_sum = re_sum + 2.0d0*(real(Xk)*cos(angle) - aimag(Xk)*sin(angle))
        end do
        k      = nm - 1
        angle  = two_pi_over_n * real(j*k, c_double)
        re_sum = re_sum + real(inp(nm, row), c_double) * cos(angle)
        out(j+1, row) = re_sum
      end do
    end do

  end subroutine irfft_rows

end module sphere_random_fields_mod
