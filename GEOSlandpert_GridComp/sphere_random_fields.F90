!
! sphere_random_fields.F90
!
! Module for generating isotropic Gaussian random fields directly on the
! sphere using the SMERFS library (Creasey & Lang 2018).
!
! Parallel strategy — Option A (m-mode distribution):
!
!   build_sphere_filter (comm optional):
!     Rank 0 computes the full all_innov / all_trans tables, then
!     MPI_Bcast to all other ranks on the node.  Every rank therefore
!     holds a copy of the complete filter.  The broadcast is cheap
!     relative to the hypergeometric evaluation cost.
!
!   draw_sphere_realisation (comm optional):
!     When comm is supplied the ranks within a node split the n_m
!     Fourier m-modes evenly (Option A).  Each rank walks the Kalman
!     filter for its strip  j_lo:j_hi  of m-modes independently (no
!     communication needed during the walk).  The resulting partial
!     half-spectrum is written into a shared MPI window
!     res_cplx(n_m, nz), then a fence ensures all ranks see the
!     complete spectrum before the IRFFT.  The IRFFT output rows are
!     also distributed: each rank computes rfield_out(:, row_lo:row_hi)
!     into the shared rfield window.  After a final fence every rank
!     returns with a complete rfield_out(nphi, nz).
!
!     Without comm the subroutine falls back to the serial path.
!
! Key public interface:
!   sphere_filter_type      - derived type holding pre-computed filter
!   build_sphere_filter     - one-time setup (rank 0 computes, bcast)
!   draw_sphere_realisation - per-call random field generation (parallel)
!   finalize_sphere_filter  - deallocate allocatables inside sf
!   sphere_to_latlon        - remap sphere field to equidistant lat-lon
!   xcorr_deg_to_c2         - derive c2 from e-folding distance in degrees
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
  use nr_ran2_gasdev, only: NRANDSEED, init_randseed, nr_ran2_2d, nr_gasdev
  use smerfs_interface
#ifdef MKL_AVAILABLE
  use MKL_DFTI
#endif

  implicit none
  private

  public :: sphere_random_fields
  public :: sphere_random_fields_id

  integer, parameter :: MORD_SPHERE = 2


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
    integer :: n_m      = 0   ! nphi/2 + 1  (half-spectrum width)
    integer :: uhalf    = 0   ! nz/2 + 1    (upper hemisphere ring count)
    integer :: Mord     = MORD_SPHERE

    ! Power-spectrum coefficients:  C_l = 1 / (c0 + c2*(l*(l+1))^2)
    real(c_double) :: c0 = 1.0d0
    real(c_double) :: c2 = 0.0d0

    ! Pre-computed state-space filter arrays (full set, all m-modes):
    !   all_innov(Mord, Mord, uhalf,   n_m)
    !   all_trans(Mord, Mord, uhalf-1, n_m)
    real(c_double), allocatable :: all_innov(:,:,:,:)
    real(c_double), allocatable :: all_trans(:,:,:,:)

    ! cos(theta) at each upper-hemisphere ring (uhalf values, equator first)
    real(c_double), allocatable :: z_pts(:)

    ! MPI context stored at build time for reuse in draw_sphere_realisation.
    ! node_comm is MPI_COMM_NULL when running serially (no comm supplied).
    integer :: node_comm   = MPI_COMM_NULL
    integer :: node_rank   = 0
    integer :: node_npes   = 1
    ! m-mode strip owned by this rank: [j_lo, j_hi]  (1-based, inclusive)
    integer :: j_lo        = 1
    integer :: j_hi        = 0   ! set to n_m in serial path

    ! Ring strip for IRFFT (set at build time, reused in draw)
    integer :: row_lo      = 1
    integer :: row_hi      = 0
    integer :: n_row_local = 0

#ifdef MKL_AVAILABLE
    ! MKL DFTI descriptor for the batched 1-D real backward FFT.
    ! One descriptor per rank, created once in build_sphere_filter.
    ! Transform: complex CCE half-spectrum (n_m) -> real (nphi),
    ! batched over n_row_local rings.
    type(DFTI_DESCRIPTOR), pointer :: dfti_irfft => null()
#endif

    logical :: initialized = .false.

  contains
  end type sphere_filter_type

  type, extends(abstract_random_fields) :: sphere_random_fields
    private
    type(sphere_filter_type) :: filter
    real(c_double), public :: c0 = 1.0d0
    real(c_double), public :: c2 = 0.0d0
    integer :: output_nlon = 0
    integer :: output_nlat = 0
    integer, public :: nz_sphere = 0
    integer, public :: nphi_sphere = 0
  contains
     procedure, public :: generate_2d_Random_field
     procedure, public :: generate_white_field
    procedure, public :: generate_2d_Random_field_seeded
    procedure, public :: sphere_to_latlon
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

  function new_sphere_random_fields(pert_param, pert_grid_f, comm, rc, c2) result(rf)
    type(sphere_random_fields) :: rf
    type(pert_param_type), intent(in) :: pert_param
    type(grid_def_type), intent(in) :: pert_grid_f
    integer, optional, intent(in) :: comm
    integer, optional, intent(out) :: rc
    real(c_double), optional, intent(in) :: c2
    real(c_double) :: xcorr_deg
    real(c_double) :: c2_use
    character(len=:), allocatable :: id_string

    xcorr_deg = real(pert_param%xcorr, c_double)
    ! Sphere perturbations always use the fine grid; no coarsening is applied.
    rf%nphi_sphere = pert_grid_f%N_lon
    rf%nz_sphere = pert_grid_f%N_lat
    rf%output_nlon = pert_grid_f%N_lon
    rf%output_nlat = pert_grid_f%N_lat
    id_string = sphere_random_fields_id(pert_param, pert_grid_f)
    rf%ID_string = id_string
    c2_use = -1.0d0
    if (present(c2)) c2_use = c2
    if (c2_use <= 0.0d0) then
      if (xcorr_deg > 0.0d0) then
        call xcorr_deg_to_c2(xcorr_deg, rf%nz_sphere, c2_use)
      else
        c2_use = 1.0d-12
      end if
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
    real :: field_mean, field_norm
    allocate(field_d(this%filter%nphi, this%filter%nz))
    allocate(field2_d(this%filter%nphi, this%filter%nz))
    call generate_2d_random_field_sphere(this%filter, rseed, field_d, field2_d, lx, ly, dx, dy)
    allocate(output_field(size(rfield,1), size(rfield,2)))
    allocate(output_field2(size(rfield2,1), size(rfield2,2)))
    call sphere_to_latlon_field(field_d, this%filter, size(rfield,1), size(rfield,2), output_field)
    call sphere_to_latlon_field(field2_d, this%filter, size(rfield2,1), size(rfield2,2), output_field2)
    rfield = output_field
    rfield2 = output_field2
    field_mean = sum(rfield) / real(size(rfield))
    rfield = rfield - field_mean
    field_norm = sqrt(sum(rfield**2) / real(size(rfield)))
    if (field_norm > 0.0) rfield = rfield / field_norm
    field_mean = sum(rfield2) / real(size(rfield2))
    rfield2 = rfield2 - field_mean
    field_norm = sqrt(sum(rfield2**2) / real(size(rfield2)))
    if (field_norm > 0.0) rfield2 = rfield2 / field_norm
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

  subroutine sphere_to_latlon(this, sphere_field, nlon_out, nlat_out, latlon_field)
    class(sphere_random_fields), intent(in) :: this
    real, intent(in) :: sphere_field(this%filter%nphi, this%filter%nz)
    integer, intent(in) :: nlon_out, nlat_out
    real, intent(out) :: latlon_field(nlon_out, nlat_out)
    real(c_double), allocatable :: sphere_d(:,:)
    real, allocatable :: latlon_d(:,:)
    allocate(sphere_d(this%filter%nphi, this%filter%nz))
    allocate(latlon_d(nlon_out, nlat_out))
    sphere_d = real(sphere_field, c_double)
    call sphere_to_latlon_field(sphere_d, this%filter, nlon_out, nlat_out, latlon_d)
    latlon_field = latlon_d
    deallocate(sphere_d, latlon_d)
  end subroutine sphere_to_latlon

  subroutine draw_sphere_realisation_seeded(sf, seed, rfield_out)
    type(sphere_filter_type), intent(in) :: sf
    integer(c_int32_t), intent(in) :: seed
    real(c_double), intent(out) :: rfield_out(sf%nphi, sf%nz)
    integer(c_int32_t), allocatable :: rand_ints(:)
    integer :: n_needed, n_got, i
    integer(c_int32_t) :: state

    n_needed = sf%nz * sf%Mord * (sf%j_hi - sf%j_lo + 1) * 2
    n_got = n_needed + n_needed/32 + 64
    allocate(rand_ints(n_got))
    state = seed
    do i = 1, n_got
      state = 1664525_c_int32_t * state + 1013904223_c_int32_t
      rand_ints(i) = state
    end do
    call draw_sphere_realisation(sf, rand_ints, n_got, rfield_out)
    deallocate(rand_ints)
  end subroutine draw_sphere_realisation_seeded

  ! ======================================================================
  !
  ! build_sphere_filter
  !
  ! Pre-compute state-space filter coefficients for the isotropic GRF
  ! with power spectrum  C_l = 1 / (c0 + c2*(l*(l+1))^2).
  !
  ! Parallel strategy:
  !   When comm is supplied, rank 0 of the global communicator runs the
  !   full computation (hyp_llp1_f loop + state_space_f loop) and then
  !   broadcasts all_innov and all_trans to all other ranks.  This is
  !   cheaper than having every rank redo the expensive hypergeometric
  !   evaluations.  The node communicator is also split here and stored
  !   in sf for reuse in draw_sphere_realisation.
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

    integer :: uhalf, n_m, m_max
    integer :: i, iroot, imode, im_hyp, rc
    integer :: global_rank, mpierr
    integer :: innov_size, trans_size

    ! Partial fraction decomposition (Mord=2 poles)
    complex(c_double_complex) :: roots_llp1(MORD_SPHERE)
    complex(c_double_complex) :: norms_pf(MORD_SPHERE)
    complex(c_double_complex) :: norm, llp1_val, lam_val, sin_pi_lam
    real(c_double)            :: llp1_re, llp1_im
    real(c_double)            :: disc_re, tmp

    ! z-points and twiddle factors
    real(c_double), allocatable :: tau_p(:,:)     ! (uhalf, 2*MORD_SPHERE-1)
    real(c_double), allocatable :: eta_ratio(:)   ! (uhalf-1)
    real(c_double), allocatable :: xvals(:)       ! (uhalf)
    real(c_double), allocatable :: yvals(:)       ! (uhalf)
    real(c_double)              :: theta_val

    ! Hypergeometric tables
    complex(c_double_complex), allocatable :: Fmat(:,:)  ! (uhalf, m_max+MORD_SPHERE+1)
    complex(c_double_complex), allocatable :: Hmat(:,:)  ! (uhalf, m_max+MORD_SPHERE+1)

    ! Covariance arrays
    real(c_double), allocatable :: cov(:,:,:,:)        ! (MORD_SPHERE,MORD_SPHERE,uhalf,0:m_max)
    real(c_double), allocatable :: cross_cov(:,:,:,:)  ! (MORD_SPHERE,MORD_SPHERE,uhalf-1,0:m_max)

    ! Temporary per-mode state-space arrays
    real(c_double), allocatable :: innov_m(:,:,:)   ! (MORD_SPHERE,MORD_SPHERE,uhalf)
    real(c_double), allocatable :: trans_m(:,:,:)   ! (MORD_SPHERE,MORD_SPHERE,uhalf-1)

    ! --------------------------------------------------------------------

    uhalf = nz/2 + 1
    n_m   = nphi/2 + 1
    m_max = n_m - 1

    ! Store grid info in sf
    sf%nz    = nz
    sf%nphi  = nphi
    sf%n_m   = n_m
    sf%uhalf = uhalf
    sf%Mord  = MORD_SPHERE
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
      ! Also get global rank to decide who does the filter computation
      call MPI_Comm_rank(comm, global_rank, mpierr)
    else
      sf%node_comm = MPI_COMM_NULL
      sf%node_rank = 0
      sf%node_npes = 1
      global_rank  = 0
    end if

    ! m-mode strip: distribute n_m modes as evenly as possible
    call compute_strip(n_m, sf%node_npes, sf%node_rank, sf%j_lo, sf%j_hi)

    ! --- Allocate filter arrays (every rank needs the full tables) ------
    allocate(sf%all_innov(MORD_SPHERE, MORD_SPHERE, uhalf,   n_m))
    allocate(sf%all_trans(MORD_SPHERE, MORD_SPHERE, uhalf-1, n_m))

    ! --- z-points -------------------------------------------------------
    allocate(sf%z_pts(uhalf))
    do i = 1, uhalf
      theta_val   = (real(uhalf - i, c_double) + 0.5d0) * &
                    (PI / real(nz, c_double))
      sf%z_pts(i) = cos(theta_val)
    end do
  
    ! --- Only rank 0 (globally) performs the expensive computation ------
    ! All other ranks wait for the MPI_Bcast below.

    if (global_rank == 0) then

      ! Twiddle factors
      allocate(tau_p   (uhalf,   2*MORD_SPHERE-1))
      allocate(eta_ratio(uhalf-1))
      tau_p(:, MORD_SPHERE)   = 1.0d0
      tau_p(:, MORD_SPHERE+1) = sqrt((1.0d0 - sf%z_pts) / (1.0d0 + sf%z_pts))
      tau_p(:, MORD_SPHERE-1) = 1.0d0 / tau_p(:, MORD_SPHERE+1)
      do i = 2, MORD_SPHERE-1
        tau_p(:, MORD_SPHERE+i) = tau_p(:, MORD_SPHERE+i-1) * tau_p(:, MORD_SPHERE+1)
        tau_p(:, MORD_SPHERE-i) = tau_p(:, MORD_SPHERE-i+1) * tau_p(:, MORD_SPHERE-1)
      end do
      do i = 1, uhalf-1
        eta_ratio(i) = tau_p(i+1, MORD_SPHERE+1) * tau_p(i, MORD_SPHERE-1)
      end do

      allocate(xvals(uhalf))
      allocate(yvals(uhalf))
      xvals = 0.5d0 * (1.0d0 - sf%z_pts)
      yvals = 0.5d0 * (1.0d0 + sf%z_pts)

      ! Partial fraction decomposition of 1/(c0 + c2*k^2), c1=0
      disc_re = -4.0d0 * c0 * c2
      if (disc_re >= 0.0d0) then
        roots_llp1(1) = cmplx( sqrt(disc_re)/(2.0d0*c2), 0.0d0, c_double_complex)
        roots_llp1(2) = cmplx(-sqrt(disc_re)/(2.0d0*c2), 0.0d0, c_double_complex)
      else
        tmp = sqrt(-disc_re) / (2.0d0*c2)
        roots_llp1(1) = cmplx(0.0d0,  tmp, c_double_complex)
        roots_llp1(2) = cmplx(0.0d0, -tmp, c_double_complex)
      end if
      do iroot = 1, MORD_SPHERE
        norms_pf(iroot) = 1.0d0 / &
          (2.0d0 * cmplx(c2, 0.0d0, c_double_complex) * roots_llp1(iroot))
      end do

      ! Covariance accumulators
      allocate(cov      (MORD_SPHERE, MORD_SPHERE, uhalf,   0:m_max))
      allocate(cross_cov(MORD_SPHERE, MORD_SPHERE, uhalf-1, 0:m_max))
      allocate(Fmat(uhalf, m_max+MORD_SPHERE+1))
      allocate(Hmat(uhalf, m_max+MORD_SPHERE+1))
      allocate(innov_m(MORD_SPHERE, MORD_SPHERE, uhalf))
      allocate(trans_m(MORD_SPHERE, MORD_SPHERE, uhalf-1))
      cov       = 0.0d0
      cross_cov = 0.0d0

      ! Hypergeometric evaluation + covariance accumulation
      do iroot = 1, MORD_SPHERE
        llp1_val = roots_llp1(iroot)
        llp1_re  = real(llp1_val, c_double)
        llp1_im  = aimag(llp1_val)
        norm     = norms_pf(iroot)

        if (llp1_re < -0.25d0) then
          lam_val = cmplx(-0.5d0, 0.0d0, c_double_complex) + &
                    cmplx(0.0d0,  1.0d0, c_double_complex) * sqrt(-0.25d0 - llp1_val)
        else
          lam_val = cmplx(-0.5d0, 0.0d0, c_double_complex) - &
                    sqrt(cmplx(0.25d0, 0.0d0, c_double_complex) + llp1_val)
        end if
        sin_pi_lam = sin(PI * lam_val)
        norm = -(0.25d0 / PI) * norm * PI / sin_pi_lam

        do im_hyp = 0, m_max + MORD_SPHERE
          call hyp_llp1_f(llp1_re, llp1_im, im_hyp, uhalf, xvals, Fmat(:,im_hyp+1), rc)
          if (rc /= 0) then
            write(*,'(a,i0)') 'build_sphere_filter: hyp_llp1_f F failed at im_hyp=', im_hyp
            stop 1
          end if
          call hyp_llp1_f(llp1_re, llp1_im, im_hyp, uhalf, yvals, Hmat(:,im_hyp+1), rc)
          if (rc /= 0) then
            write(*,'(a,i0)') 'build_sphere_filter: hyp_llp1_f H failed at im_hyp=', im_hyp
            stop 1
          end if
        end do

        call update_cov_f(m_max, uhalf, MORD_SPHERE,        &
                          real(norm,c_double), aimag(norm),  &
                          llp1_re, llp1_im,                  &
                          Fmat, Hmat, tau_p, eta_ratio,      &
                          cov, cross_cov, rc)
        if (rc /= 0) then
          write(*,'(a,i0)') 'build_sphere_filter: update_cov_f failed, rc=', rc
          stop 1
        end if
      end do  ! iroot

      ! State-space decomposition for each m-mode
      do imode = 0, n_m-1
        call state_space_f(uhalf, MORD_SPHERE,              &
                           cross_cov(:,:,:,imode),          &
                           cov(:,:,:,imode),                &
                           innov_m, trans_m, rc)
        if (rc /= 0) then
          write(*,'(a,i0,a,i0)') &
            'build_sphere_filter: state_space_f failed at imode=', imode, ', rc=', rc
          stop 1
        end if
        sf%all_innov(:,:,:,imode+1) = innov_m
        sf%all_trans(:,:,:,imode+1) = trans_m
      end do

      deallocate(tau_p, eta_ratio, xvals, yvals)
      deallocate(Fmat, Hmat, cov, cross_cov, innov_m, trans_m)

    end if  ! global_rank == 0

    ! --- Broadcast filter tables to all ranks ---------------------------
    ! Use the global communicator (comm if present, else no-op since
    ! global_rank == 0 is the only rank).
    if (present(comm)) then
      innov_size = MORD_SPHERE * MORD_SPHERE * uhalf   * n_m
      trans_size = MORD_SPHERE * MORD_SPHERE * (uhalf-1) * n_m
      call MPI_Bcast(sf%all_innov, innov_size, MPI_DOUBLE_PRECISION, 0, comm, mpierr)
      call MPI_Bcast(sf%all_trans, trans_size, MPI_DOUBLE_PRECISION, 0, comm, mpierr)
      call MPI_Bcast(sf%z_pts,     uhalf,      MPI_DOUBLE_PRECISION, 0, comm, mpierr)
    end if

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
    ! Storage:   DFTI_NOT_INPLACE, separate input/output buffers.
    !            Input  = complex(c_float) half-spectrum, length n_m per ring.
    !            Output = real(c_float),                  length nphi per ring.
    !
    ! Batching over n_row_local rings:
    !   DFTI_NUMBER_OF_TRANSFORMS = n_row_local
    !   DFTI_INPUT_DISTANCE  = n_m   (complex elements between ring starts)
    !   DFTI_OUTPUT_DISTANCE = nphi  (real    elements between ring starts)
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
        ! In-place real backward FFT: input and output share the same buffer.
        ! work(:) is real(c_float), size 2*n_m*n_row_local.
        ! MKL CCE in-place format:
        !   INPUT_DISTANCE  = distance between ring starts in complex elements = n_m
        !   OUTPUT_DISTANCE = distance between ring starts in real elements    = 2*n_m
        dfti_stat = DftiSetValue(sf%dfti_irfft, &
                                 DFTI_NUMBER_OF_TRANSFORMS, sf%n_row_local)
        dfti_stat = DftiSetValue(sf%dfti_irfft, DFTI_INPUT_DISTANCE,  n_m)
        dfti_stat = DftiSetValue(sf%dfti_irfft, DFTI_OUTPUT_DISTANCE, 2*n_m)
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

  end subroutine build_sphere_filter

  ! ======================================================================
  !
  ! draw_sphere_realisation  (Option A parallel implementation)
  !
  ! Generate one realisation of the isotropic GRF on the sphere.
  !
  ! Parallel strategy (Option A — m-mode distribution):
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
  !   subroutine runs the original serial path without any MPI calls.
  !
  ! Arguments:
  !   sf           (in)  : initialised sphere_filter_type
  !   rand_ints    (in)  : array of random 32-bit integers
  !                        Serial:   size >= n_rand_got (all modes)
  !                        Parallel: size >= strip_rand_got (local strip)
  !   n_rand_got   (in)  : length of rand_ints supplied
  !   rfield_out   (out) : real(c_double) field (sf%nphi, sf%nz)
  !                        On return every rank holds the complete field.

  subroutine draw_sphere_realisation(sf, rand_ints, n_rand_got, rfield_out)

    type(sphere_filter_type),            intent(in)  :: sf
    integer,                             intent(in)  :: n_rand_got
    integer(c_int32_t),                  intent(in)  :: rand_ints(n_rand_got)
    real(c_double),                      intent(out) :: rfield_out(sf%nphi, sf%nz)

    ! --- local variables ------------------------------------------------

    integer :: uhalf, n_m, nz, nphi, Mord
    integer :: j_lo, j_hi, n_m_local          ! m-mode strip for this rank
    integer :: row_lo, row_hi, n_row_local     ! ring strip for IRFFT (from sf)
    integer :: n_rand_needed_local
     integer :: i, j, ip, iq, nn, rc
     integer :: skip
    logical :: parallel

    real(c_float), allocatable :: noise_real(:)          ! standard normals
    complex(c_float), allocatable :: noise(:,:,:)        ! (Mord, n_m_local, nz)
     complex(c_float), allocatable :: fp_f(:,:)           ! (Mord, n_m_local)
     complex(c_float), allocatable :: fp_fstart(:,:)      ! (Mord, n_m_local)
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

    ! --------------------------------------------------------------------

    uhalf = sf%uhalf
    n_m   = sf%n_m
    nz    = sf%nz
    nphi  = sf%nphi
    Mord  = sf%Mord
    j_lo  = sf%j_lo
    j_hi  = sf%j_hi

    parallel = (sf%node_comm /= MPI_COMM_NULL)

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
      ! Serial: use private allocatables and point to them
      allocate(res_cplx_priv(n_m, nz))
      allocate(rfield_priv  (nphi, nz))
      res_cplx => res_cplx_priv
      rfield_win => rfield_priv
      n_m_local = n_m    ! serial owns all modes
      j_lo = 1 ; j_hi = n_m
      row_lo = 1 ; row_hi = nz ; n_row_local = nz
    end if

    ! --- 2. Noise generation for the local m-mode strip -----------------
    ! n_rand_needed for this rank: nz * Mord * n_m_local * 2
    ! The south-hemisphere walk reuses noise at indices uhalf+1..nz, so
    ! the full nz rings of noise are required (not just uhalf).
    n_rand_needed_local = nz * Mord * n_m_local * 2

    allocate(noise_real(n_rand_needed_local))
    allocate(noise    (Mord, n_m_local, nz))
     allocate(fp_f     (Mord, n_m_local))
     allocate(fp_fstart(Mord, n_m_local))
     allocate(fp_fprev (Mord, n_m_local))

    if (n_m_local > 0) then
      call zigg_f(n_rand_needed_local, n_rand_got, rand_ints, noise_real, rc)
      if (rc /= 0) then
        write(*,'(a,i0,a)') &
          'WARNING draw_sphere_realisation: zigg exhausted, ', rc, ' samples missing'
      end if

      ! Pack into complex noise: noise(ip, jlocal, i) = (re + i*im) / sqrt(2)
      ! Loop over all nz rings so south-hemisphere indices (uhalf+1..nz)
      ! are filled from the same noise stream as the reference Python code.
      nn = 0
      do i = 1, nz
        do ip = 1, Mord
          do j = 1, n_m_local
            nn = nn + 1
            noise(ip, j, i) = cmplx(noise_real(2*nn-1), noise_real(2*nn), c_float) &
                              * real(0.5d0**0.5d0, c_float)
          end do
        end do
      end do

      ! --- 3. Kalman walk for local m-modes j_lo:j_hi -------------------
      ! Upward walk: equator (ring 1 in sf%z_pts) -> pole (ring uhalf)

      ! Initialise at ring 1
      do j = 1, n_m_local
        do ip = 1, Mord
          fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
          do iq = 1, Mord
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
          do ip = 1, Mord
            fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
            do iq = 1, Mord
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
      do ip = 1, Mord
        if (mod(ip-1, 2) == 0) then
             fp_f(ip, :) = fp_fstart(ip, :)
        else
          fp_f(ip, :) = -fp_fstart(ip, :)
        end if
      end do
       fp_fprev = fp_f

       ! For even nz, the equatorial transition was already represented by
       ! the northern walk, so skip that transition on the southward walk.
       skip = 1 - mod(nz, 2)

      do i = 0, nz - uhalf - 1
        do j = 1, n_m_local
          do ip = 1, Mord
            fp_f(ip, j) = cmplx(0.0, 0.0, c_float)
            do iq = 1, Mord
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

    ! --- 4. Fence: all ranks have written their m-mode strip ------------
    if (parallel) then
      call MPI_Win_fence(0, win_res, mpierr)
    end if

    ! --- 5. Distributed IRFFT: each rank transforms its ring strip ------
    ! MKL path: batched 1-D real backward FFT using the pre-built descriptor.
    ! Fallback: pure-Fortran O(nphi^2) DFT (NR path or n_row_local==0).
#ifdef MKL_AVAILABLE
    if (n_row_local > 0) then
      block
        ! In-place real backward FFT using the pre-built descriptor.
        ! The work buffer holds n_row_local rings, each padded to 2*n_m
        ! real(c_float) elements (= n_m complex elements).
        ! MKL reads the complex half-spectrum from the first n_m complex
        ! positions (= 2*n_m reals) and writes the real output (nphi reals)
        ! starting at the same address.
        real(c_float), allocatable, target :: work(:)   ! (2*n_m * n_row_local)
        type(c_ptr) :: cptr_w
        complex(c_float), pointer :: cplx_view(:,:)  ! (n_m, n_row_local)
        real(c_float),    pointer :: real_view(:,:)   ! (2*n_m, n_row_local)
        integer :: dfti_stat, row

        allocate(work(2 * n_m * n_row_local))

        ! Map complex view onto work buffer, copy half-spectrum in
        cptr_w = c_loc(work(1))
        call c_f_pointer(cptr_w, cplx_view, [n_m, n_row_local])
        cplx_view = res_cplx(:, row_lo:row_hi)

        ! In-place backward FFT
        dfti_stat = DftiComputeBackward(sf%dfti_irfft, work)
        if (dfti_stat /= DFTI_NO_ERROR) then
          write(*,'(a,i0)') &
            'draw_sphere_realisation: DftiComputeBackward failed, stat=', dfti_stat
          stop 1
        end if

        ! Real output occupies first nphi reals of each padded row (stride 2*n_m)
        call c_f_pointer(cptr_w, real_view, [2*n_m, n_row_local])
        do row = 1, n_row_local
          rfield_win(:, row_lo + row - 1) = real(real_view(1:nphi, row), c_double)
        end do

        deallocate(work)
      end block
    end if
#else
    if (n_row_local > 0) then
      call irfft_rows(res_cplx(:, row_lo:row_hi), n_row_local, n_m, nphi, &
                      rfield_win(:, row_lo:row_hi))
    end if
#endif

    ! --- 6. Fence: all ranks have written their ring strip ---------------
    if (parallel) call MPI_Win_fence(0, win_rf, mpierr)

    ! --- 7. Normalise to unit variance (global empirical) ----------------
    ! The IRFFT output has a variance that depends on nz, nphi, Mord, and
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

  end subroutine draw_sphere_realisation

  ! ======================================================================
  !
  ! finalize_sphere_filter

  subroutine finalize_sphere_filter(sf)
    type(sphere_filter_type), intent(inout) :: sf
    integer :: mpierr

    if (allocated(sf%all_innov)) deallocate(sf%all_innov)
    if (allocated(sf%all_trans)) deallocate(sf%all_trans)
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
  ! sphere_to_latlon
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
  ! xcorr_deg is interpreted as the Gaussian sigma of the covariance, i.e.
  ! the angle at which C(theta) = C(0) * exp(-0.5).  This matches the
  ! flat-FFT path in random_fields.F90, where the power spectrum amplitude
  ! is exp(-0.25*(lx*kx)^2), giving covariance C(r) = exp(-r^2/(2*lx^2))
  ! and Gaussian sigma = lx = xcorr [in whatever units lx and dx share].
  !
  ! The true e-folding distance (where C = C(0)/e) is xcorr_deg * sqrt(2).
  !
  ! Uses bisection on the analytic Legendre-series covariance.

  subroutine xcorr_deg_to_c2(xcorr_deg, nz, c2_out)

    real(c_double), intent(in)  :: xcorr_deg
    integer,        intent(in)  :: nz
    real(c_double), intent(out) :: c2_out

    integer,        parameter :: lmax_mult = 3
    real(c_double), parameter :: tol_deg = 0.05d0
    integer,        parameter :: max_iter = 80

    integer        :: iter, lmax
    real(c_double) :: c2_lo, c2_hi, c2_mid, theta_mid, xcorr_rad

    xcorr_rad = xcorr_deg * MAPL_DEGREES_TO_RADIANS_R8
    lmax      = lmax_mult * nz

    ! Bracket: larger c2 -> smoother spectrum (more low-l power) -> LONGER
    ! correlation.  So c2_lo gives SHORT correlation, c2_hi gives LONG.
    ! c2_lo = 1e-10 -> Gaussian sigma ~ 3.5 deg  (very short)
    ! c2_hi = 1e+2  -> Gaussian sigma ~ 127 deg  (nearly uniform)
    c2_lo = 1.0d-10
    c2_hi = 1.0d2

    if (xcorr_rad < compute_efold_theta(c2_lo, lmax)) then
      write(*,'(a,f8.3,a)') 'xcorr_deg_to_c2: xcorr=', xcorr_deg, &
        ' deg is very small (< 5 deg); using c2_lo'
      c2_out = c2_lo ; return
    end if
    if (xcorr_rad > compute_efold_theta(c2_hi, lmax)) then
      write(*,'(a,f8.3,a)') 'xcorr_deg_to_c2: xcorr=', xcorr_deg, &
        ' deg is very large; using c2_hi'
      c2_out = c2_hi ; return
    end if

    ! Bisection: theta_mid > xcorr_rad means c2_mid gives too long a scale,
    ! so move the upper bracket down.
    c2_mid = c2_lo
    do iter = 1, max_iter
      c2_mid    = 0.5d0 * (c2_lo + c2_hi)
      theta_mid = compute_efold_theta(c2_mid, lmax)
      if (abs(theta_mid - xcorr_rad) < tol_deg * MAPL_DEGREES_TO_RADIANS_R8) exit
      if (theta_mid > xcorr_rad) then
        c2_hi = c2_mid   ! c2_mid gives too long a scale; reduce c2
      else
        c2_lo = c2_mid   ! c2_mid gives too short a scale; increase c2
      end if
    end do
    c2_out = c2_mid

  contains

    real(c_double) function compute_efold_theta(c2_val, lmax_in)
      real(c_double), intent(in) :: c2_val
      integer,        intent(in) :: lmax_in

      real(c_double) :: C0_z1, C_z, z, dz, C_l_val, llp1, P0, P1, P2
      integer        :: l, iz, nz_scan
      real(c_double), parameter :: inv4pi = 0.25d0 / PI

      nz_scan = 500
      dz = 2.0d0 / real(nz_scan, c_double)

      C0_z1 = 0.0d0 ; P0 = 1.0d0 ; P1 = 1.0d0
      do l = 0, lmax_in
        llp1    = real(l,c_double) * real(l+1,c_double)
        C_l_val = 1.0d0 / (1.0d0 + c2_val * llp1**2)
        if (l == 0) then
          C0_z1 = C0_z1 + C_l_val * inv4pi
        else if (l == 1) then
          C0_z1 = C0_z1 + C_l_val * 3.0d0 * inv4pi
        else
          P2 = (real(2*l-1,c_double)*P1 - real(l-1,c_double)*P0) / real(l,c_double)
          C0_z1 = C0_z1 + C_l_val * real(2*l+1,c_double) * inv4pi
          P0 = P1 ; P1 = P2
        end if
      end do

      compute_efold_theta = PI
      do iz = 0, nz_scan
        z = 1.0d0 - real(iz,c_double) * dz
        C_z = 0.0d0 ; P0 = 1.0d0 ; P1 = z
        do l = 0, lmax_in
          llp1    = real(l,c_double) * real(l+1,c_double)
          C_l_val = 1.0d0 / (1.0d0 + c2_val * llp1**2)
          if (l == 0) then
            C_z = C_z + C_l_val * inv4pi
          else if (l == 1) then
            C_z = C_z + C_l_val * 3.0d0 * z * inv4pi
          else
            P2 = (real(2*l-1,c_double)*z*P1 - real(l-1,c_double)*P0) / real(l,c_double)
            C_z = C_z + C_l_val * real(2*l+1,c_double) * P2 * inv4pi
            P0 = P1 ; P1 = P2
          end if
        end do
        ! Use exp(-0.5) as the target so that xcorr_deg is the Gaussian sigma
        ! (C(sigma) = C(0)*exp(-0.5)), matching the flat-FFT definition where
        ! the covariance is exp(-r^2/(2*xcorr^2)).
        if (C_z <= exp(-0.5d0) * C0_z1) then
          compute_efold_theta = acos(max(-1.0d0, min(1.0d0, z)))
          return
        end if
      end do
    end function compute_efold_theta

  end subroutine xcorr_deg_to_c2

  subroutine generate_2d_random_field_sphere(this, rseed, rfield, rfield2, lx, ly, dx, dy)
    class(sphere_filter_type), intent(inout) :: this
    integer, intent(inout) :: rseed(NRANDSEED)
    real(c_double), intent(out) :: rfield(this%nphi, this%nz)
    real(c_double), intent(out) :: rfield2(this%nphi, this%nz)
    real, intent(in) :: lx, ly, dx, dy
    integer(c_int32_t), allocatable :: ints(:)
    real, allocatable :: uniforms(:)
    integer :: n_needed, n_got

    n_needed = this%nz * this%Mord * (this%j_hi - this%j_lo + 1) * 2
    n_got = n_needed + n_needed/32 + 64
    allocate(ints(n_got), uniforms(n_got))
    call nr_ran2_2d(1, n_got, rseed, uniforms)
    ints = int(uniforms * 2147483647.0, c_int32_t)
    call draw_sphere_realisation(this, ints, n_got, rfield)
    call nr_ran2_2d(1, n_got, rseed, uniforms)
    ints = int(uniforms * 2147483647.0, c_int32_t)
    call draw_sphere_realisation(this, ints, n_got, rfield2)
    deallocate(ints, uniforms)
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
    real :: field_mean, field_norm

    n_needed = this%filter%nz * this%filter%Mord * &
               (this%filter%j_hi - this%filter%j_lo + 1) * 2
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
    field_mean = sum(rfield) / real(size(rfield))
    rfield = rfield - field_mean
    field_norm = sqrt(sum(rfield**2) / real(size(rfield)))
    if (field_norm > 0.0) rfield = rfield / field_norm
    field_mean = sum(rfield2) / real(size(rfield2))
    rfield2 = rfield2 - field_mean
    field_norm = sqrt(sum(rfield2**2) / real(size(rfield2)))
    if (field_norm > 0.0) rfield2 = rfield2 / field_norm
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
  ! Private: analytic_var
  !
  ! Compute the analytic lag-0 variance of the GRF:
  !   C(0) = sum_{l=0}^{lmax} (2l+1) * C_l / (4*pi)
  !   C_l = 1 / (c0 + c2*(l*(l+1))^2)

  real(c_double) function analytic_var(c0, c2, lmax)
    real(c_double), intent(in) :: c0, c2
    integer,        intent(in) :: lmax
    integer        :: l
    real(c_double) :: llp1, Cl
    real(c_double), parameter :: inv4pi = 0.25d0 / PI
    analytic_var = 0.0d0
    do l = 0, lmax
      llp1 = real(l,c_double) * real(l+1,c_double)
      Cl   = 1.0d0 / (c0 + c2 * llp1**2)
      analytic_var = analytic_var + real(2*l+1, c_double) * Cl * inv4pi
    end do
  end function analytic_var

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

  ! ======================================================================
  ! Private: native ring z coordinate

  real(c_double) function ring_z(k, sf)
    integer, intent(in) :: k
    type(sphere_filter_type), intent(in) :: sf

    if (k <= sf%uhalf) then
      ring_z = sf%z_pts(sf%uhalf + 1 - k)
    else
      ring_z = -sf%z_pts(k - sf%uhalf + 1)
    end if
  end function ring_z

  ! Private: find_nearest_ring

  integer function find_nearest_ring(z_target, sf)
    real(c_double),           intent(in) :: z_target
    type(sphere_filter_type), intent(in) :: sf

    integer        :: k, k_best
    real(c_double) :: dz, dz_best, z_k

    ! Ring k layout:
    !   k = 1..uhalf   : north walk result, z = +z_pts(uhalf+1-k)  [pole->equator]
    !   k = uhalf+1..nz: south walk result, z = -z_pts(k-uhalf+1)  [equator->S-pole]
    !
    ! The south-walk step i_s (0-based) stores in ring k = uhalf+i_s+1 the field
    ! produced by the filter transition from z_pts(i_s+skip+1) to z_pts(i_s+skip+2),
    ! i.e. the field lives at z = -z_pts(i_s+skip+2) = -z_pts(k-uhalf+1)  (skip=1 for
    ! even nz).  The previous formula used -z_pts(k-uhalf), which was off by one.
    dz_best = 1.0d10 ; k_best = 1
    do k = 1, sf%nz
      z_k = ring_z(k, sf)
      dz = abs(z_k - z_target)
      if (dz < dz_best) then
        dz_best = dz ; k_best = k
      end if
    end do
    find_nearest_ring = k_best
  end function find_nearest_ring

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
