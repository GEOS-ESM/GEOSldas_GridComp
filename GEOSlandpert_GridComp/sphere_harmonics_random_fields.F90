!
! sphere_harmonics_random_fields.F90
!
! Isotropic Gaussian random fields on the sphere by direct spherical-
! harmonic (SH) synthesis, with power spectrum
!
!     C_l = 1 / (c0 + c2 * (l(l+1))^2),   c0 = 1,
!
! c2 calibrated from xcorr (Gaussian sigma, degrees) by xcorr_deg_to_c2.
!
! Synthesis on the regular lat-lon grid (ring centres theta_j=(j-0.5)*pi/nlat):
!
!   f(theta,phi) = sum_l [ a_l0 Y_l0 + 2 Re sum_{m>0} a_lm Y_lm ]
!   a_l0 ~ N(0,C_l),  Re/Im a_lm ~ N(0,C_l/2),  l,m <= L = nlat-1
!
! Build (constructor): c2, sqrt(C_l), and (optionally) the table of
!   normalised associated Legendre functions Pbar_lm(cos theta_j) for the
!   north-hemisphere rings and the m's owned by this rank (snake order).
! Draw: random a_lm -> per-m DGEMM against Pbar (even/odd l split gives both
!   hemispheres) -> MPI_Reduce_scatter to ring blocks -> batched inverse
!   real FFT -> MPI_Allgatherv.  Two fields are drawn per call.
!
! The random coefficients are taken from one global stream (identical on
! every rank) so the result is independent of the number of ranks and the
! rseed is advanced identically everywhere.
!
! Public interface:
!   sphere_harmonics_random_fields    - abstract_random_fields implementation
!   sphere_harmonics_random_fields_id - cache key for the land_pert generator map
!
! -----------------------------------------------------------------------

module sphere_harmonics_random_fields_mod
  use abstract_random_fieldsMod, only: abstract_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  use MAPL_MathConstantsMod, only: PI => MAPL_PI_R8, MAPL_DEGREES_TO_RADIANS_R8
  use, intrinsic :: iso_fortran_env, only: REAL64, INT64
  use mpi
  use nr_ran2_gasdev, only: NRANDSEED, nr_ran2_2d, nr_gasdev, init_randseed
#ifdef MKL_AVAILABLE
  use MKL_DFTI
#endif

  implicit none
  private

  public :: sphere_harmonics_random_fields
  public :: sphere_harmonics_random_fields_id

  type, extends(abstract_random_fields) :: sphere_harmonics_random_fields
    private
    ! MPI context (comm == MPI_COMM_NULL -> serial)
    integer :: comm = MPI_COMM_NULL
    integer :: rank = 0, npes = 1
    ! Grid / truncation
    integer :: L = 0, nh = 0, nm = 0
    logical :: use_table = .false.
    real(REAL64) :: scale = 1.0d0
    ! m-modes owned by this rank and their offsets
    integer :: nmloc = 0
    integer,    allocatable :: mloc(:)
    integer(INT64), allocatable :: off(:)        ! offsets into ptab
    integer(INT64), allocatable :: coff(:)       ! global coefficient offset per m (0:L)
    integer(INT64) :: ncoef = 0                  ! total (l,m) pairs, m<=l<=L
    real(REAL64), allocatable :: ptab(:)
    real(REAL64), allocatable :: sqrtC(:), lfac(:), x(:), logs(:)
    ! Ring distribution
    integer, allocatable :: rcount(:), rdispl(:)
    integer :: nloc = 0
#ifdef MKL_AVAILABLE
    type(DFTI_DESCRIPTOR), pointer :: desc => null()
#endif
    logical :: is_initialized = .false.

    real(REAL64), public :: c0 = 1.0d0
    real(REAL64), public :: c2 = 0.0d0
    integer, public :: nz_sphere = 0
    integer, public :: nphi_sphere = 0
    integer(INT64), public :: table_bytes = 0    ! global Legendre table size
  contains
    procedure, public :: generate_2d_Random_field
    procedure, public :: generate_white_field
    procedure, public :: generate_2d_Random_field_seeded
    procedure, public :: finalize
    procedure, public :: initialized
  end type sphere_harmonics_random_fields

  interface sphere_harmonics_random_fields
    module procedure new_sphere_harmonics_random_fields
  end interface sphere_harmonics_random_fields

contains

  function sphere_harmonics_random_fields_id(pert_param, pert_grid_f) result(id_string)
    type(pert_param_type), intent(in) :: pert_param
    type(grid_def_type), intent(in) :: pert_grid_f
    character(len=:), allocatable :: id_string
    character(len=256) :: id

    write(id,'(a,":",i0,":",i0,":",es16.8)') 'sphere_sh', &
         pert_grid_f%N_lat, pert_grid_f%N_lon, pert_param%xcorr
    id_string = trim(id)
  end function sphere_harmonics_random_fields_id

  ! ----------------------------------------------------------------------
  ! Constructor.  The optional
  ! use_table (default .false.) selects the precomputed Legendre table
  ! (memory ~ 2*nlat^3 bytes in total) instead of on-the-fly recursion.
  ! The table is not faster: the 4-column DGEMM is memory-bound, so
  ! streaming the table costs more than recomputing Pbar in cache.

  function new_sphere_harmonics_random_fields(pert_param, pert_grid_f, comm, rc, use_table) result(rf)
    type(sphere_harmonics_random_fields) :: rf
    type(pert_param_type), intent(in) :: pert_param
    type(grid_def_type), intent(in) :: pert_grid_f
    integer, optional, intent(in) :: comm
    integer, optional, intent(out) :: rc
    logical, optional, intent(in) :: use_table

    integer :: ierr, i, m, r, l, nl, nlat, nlon
    integer(INT64) :: tot
    real(REAL64) :: xcorr_deg, var, th, llp1, c2_use

    if (present(rc)) rc = 0
    nlon = pert_grid_f%N_lon
    nlat = pert_grid_f%N_lat
    rf%nphi_sphere = nlon
    rf%nz_sphere   = nlat
    rf%ID_string   = sphere_harmonics_random_fields_id(pert_param, pert_grid_f)
    rf%use_table   = .false.
    if (present(use_table)) rf%use_table = use_table

    if (present(comm)) then
      call MPI_Comm_dup(comm, rf%comm, ierr)
      call MPI_Comm_rank(rf%comm, rf%rank, ierr)
      call MPI_Comm_size(rf%comm, rf%npes, ierr)
    else
      rf%comm = MPI_COMM_NULL
      rf%rank = 0
      rf%npes = 1
    end if

    rf%L  = nlat - 1
    rf%nh = nlat / 2
    rf%nm = nlon / 2 + 1
    if (mod(nlat, 2) /= 0 .or. mod(nlon, 2) /= 0 .or. rf%nm < rf%L + 1) then
      if (rf%rank == 0) write(*,'(a,i0,a,i0,a)') &
        'sphere_harmonics_random_fields: need even nlat, nlon >= 2*nlat-2 (got ', &
        nlon, 'x', nlat, ')'
      if (present(rc)) then
        rc = 1; return
      end if
      stop 1
    end if

    ! --- spectrum -------------------------------------------------------
    xcorr_deg = real(pert_param%xcorr, REAL64)
    if (rf%rank == 0) then
      if (xcorr_deg > 0.0d0) then
        call xcorr_deg_to_c2(xcorr_deg, nlat, c2_use)
      else
        c2_use = 1.0d-12
      end if
    end if
    if (rf%comm /= MPI_COMM_NULL) &
      call MPI_Bcast(c2_use, 1, MPI_DOUBLE_PRECISION, 0, rf%comm, ierr)
    rf%c2 = c2_use

    allocate(rf%sqrtC(0:rf%L))
    var = 0.0d0
    do l = 0, rf%L
      llp1 = real(l,REAL64)*real(l+1,REAL64)
      rf%sqrtC(l) = sqrt(1.0d0 / (rf%c0 + rf%c2*llp1**2))
      var = var + real(2*l+1,REAL64) * rf%sqrtC(l)**2 / (4.0d0*PI)
    end do
    rf%scale = 1.0d0 / sqrt(var)

    ! --- log of Pbar_mm prefactor: log(sqrt(1/4pi) prod sqrt((2k+1)/2k)) --
    allocate(rf%lfac(0:rf%L))
    rf%lfac(0) = 0.5d0*log(1.0d0/(4.0d0*PI))
    do m = 1, rf%L
      rf%lfac(m) = rf%lfac(m-1) + 0.5d0*log(real(2*m+1,REAL64)/real(2*m,REAL64))
    end do
    allocate(rf%x(rf%nh), rf%logs(rf%nh))
    do i = 1, rf%nh
      th = (real(i,REAL64) - 0.5d0) * PI / real(nlat,REAL64)
      rf%x(i) = cos(th); rf%logs(i) = log(sin(th))
    end do

    ! --- global coefficient offsets (rank-independent random stream) ----
    allocate(rf%coff(0:rf%L))
    rf%coff(0) = 0
    do m = 1, rf%L
      rf%coff(m) = rf%coff(m-1) + int(rf%L - m + 2, INT64)
    end do
    rf%ncoef = rf%coff(rf%L) + 1

    ! --- m distribution: snake order for load balance -------------------
    allocate(rf%mloc(rf%L/rf%npes + 2))
    rf%nmloc = 0
    do m = 0, rf%L
      r = mod(m, 2*rf%npes)
      if (r >= rf%npes) r = 2*rf%npes - 1 - r
      if (r == rf%rank) then
        rf%nmloc = rf%nmloc + 1
        rf%mloc(rf%nmloc) = m
      end if
    end do

    ! --- Legendre table -------------------------------------------------
    allocate(rf%off(rf%nmloc + 1))
    rf%off(1) = 1
    do i = 1, rf%nmloc
      nl = rf%L - rf%mloc(i) + 1
      rf%off(i+1) = rf%off(i) + int(rf%nh,INT64)*int(nl,INT64)
    end do
    rf%table_bytes = 0
    if (rf%use_table) then
      allocate(rf%ptab(rf%off(rf%nmloc+1) - 1))
      do i = 1, rf%nmloc
        nl = rf%L - rf%mloc(i) + 1
        call plm_block(rf, rf%mloc(i), nl, rf%ptab(rf%off(i):rf%off(i+1)-1))
      end do
      rf%table_bytes = 8_INT64 * (rf%off(rf%nmloc+1) - 1)
    end if
    tot = rf%table_bytes
    if (rf%comm /= MPI_COMM_NULL) &
      call MPI_Allreduce(MPI_IN_PLACE, tot, 1, MPI_INTEGER8, MPI_SUM, rf%comm, ierr)
    rf%table_bytes = tot

    ! --- ring distribution ---------------------------------------------
    allocate(rf%rcount(rf%npes), rf%rdispl(rf%npes))
    do r = 0, rf%npes-1
      rf%rcount(r+1) = nlat/rf%npes + merge(1, 0, r < mod(nlat, rf%npes))
    end do
    rf%rdispl(1) = 0
    do r = 2, rf%npes
      rf%rdispl(r) = rf%rdispl(r-1) + rf%rcount(r-1)
    end do
    rf%nloc = rf%rcount(rf%rank+1)

#ifdef MKL_AVAILABLE
    block
      integer :: st
      st = DftiCreateDescriptor(rf%desc, DFTI_DOUBLE, DFTI_REAL, 1, nlon)
      st = DftiSetValue(rf%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      st = DftiSetValue(rf%desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
      st = DftiSetValue(rf%desc, DFTI_NUMBER_OF_TRANSFORMS, max(1, 2*rf%nloc))
      st = DftiSetValue(rf%desc, DFTI_INPUT_DISTANCE,  rf%nm)
      st = DftiSetValue(rf%desc, DFTI_OUTPUT_DISTANCE, nlon)
      st = DftiCommitDescriptor(rf%desc)
      if (st /= DFTI_NO_ERROR) then
        write(*,'(a,i0)') 'sphere_harmonics_random_fields: DFTI commit failed ', st
        if (present(rc)) then
          rc = st; return
        end if
        stop 1
      end if
    end block
#endif

    rf%is_initialized = .true.
  end function new_sphere_harmonics_random_fields

  ! ----------------------------------------------------------------------

  subroutine generate_2d_Random_field(this, rseed, rfield, rfield2, lx, ly, dx, dy)
    class(sphere_harmonics_random_fields), intent(inout) :: this
    integer, intent(inout) :: rseed(NRANDSEED)
    real, intent(out) :: rfield(:,:), rfield2(:,:)
    real, intent(in) :: lx, ly, dx, dy
    real, allocatable :: f1(:,:), f2(:,:)

    allocate(f1(this%nphi_sphere, this%nz_sphere), f2(this%nphi_sphere, this%nz_sphere))
    call sh_draw(this, rseed, f1, f2)
    call remap_latlon(f1, rfield)
    call remap_latlon(f2, rfield2)
    call normalize_unit_variance(rfield)
    call normalize_unit_variance(rfield2)
  end subroutine generate_2d_Random_field

  subroutine generate_white_field(this, rseed, rfield)
    class(sphere_harmonics_random_fields), intent(inout) :: this
    integer, intent(inout) :: rseed(NRANDSEED)
    real, target, intent(out) :: rfield(:,:)
    real :: values(2)
    integer :: i, j

    do j = 1, size(rfield, 2)
      do i = 1, size(rfield, 1), 2
        call nr_gasdev(rseed, values)
        rfield(i,j) = values(1)
        if (i + 1 <= size(rfield, 1)) rfield(i + 1,j) = values(2)
      end do
    end do
  end subroutine generate_white_field

  ! Integer-seeded variant: initialise a private ran2 stream from seed.
  subroutine generate_2d_Random_field_seeded(this, seed, rfield, rfield2, lx, ly, dx, dy)
    class(sphere_harmonics_random_fields), intent(inout) :: this
    integer, intent(in) :: seed
    real, intent(out) :: rfield(:,:), rfield2(:,:)
    real, intent(in) :: lx, ly, dx, dy
    integer :: rseed(NRANDSEED)

    rseed = 0
    rseed(1) = -max(abs(seed), 1)
    call init_randseed(rseed)
    call this%generate_2d_Random_field(rseed, rfield, rfield2, lx, ly, dx, dy)
  end subroutine generate_2d_Random_field_seeded

  subroutine finalize(this, rc)
    class(sphere_harmonics_random_fields), intent(inout) :: this
    integer, optional, intent(out) :: rc
    integer :: ierr
#ifdef MKL_AVAILABLE
    if (associated(this%desc)) then
      ierr = DftiFreeDescriptor(this%desc)
      nullify(this%desc)
    end if
#endif
    if (allocated(this%ptab))   deallocate(this%ptab)
    if (allocated(this%mloc))   deallocate(this%mloc)
    if (allocated(this%off))    deallocate(this%off)
    if (allocated(this%coff))   deallocate(this%coff)
    if (allocated(this%sqrtC))  deallocate(this%sqrtC)
    if (allocated(this%lfac))   deallocate(this%lfac)
    if (allocated(this%x))      deallocate(this%x)
    if (allocated(this%logs))   deallocate(this%logs)
    if (allocated(this%rcount)) deallocate(this%rcount)
    if (allocated(this%rdispl)) deallocate(this%rdispl)
    if (this%comm /= MPI_COMM_NULL) then
      call MPI_Comm_free(this%comm, ierr)
      this%comm = MPI_COMM_NULL
    end if
    this%is_initialized = .false.
    if (present(rc)) rc = 0
  end subroutine finalize

  logical function initialized(this)
    class(sphere_harmonics_random_fields), intent(in) :: this
    initialized = this%is_initialized
  end function initialized

  ! ======================================================================
  ! Pbar_lm(x_j) for l = m..L, north rings j = 1..nh, stored P(nh, nl).

  subroutine plm_block(sh, m, nl, P)
    type(sphere_harmonics_random_fields), intent(in) :: sh
    integer, intent(in) :: m, nl
    real(REAL64), intent(out) :: P(sh%nh, nl)
    integer :: j, l, k
    real(REAL64) :: a, b, rl, rm

    rm = real(m,REAL64)
    do j = 1, sh%nh
      P(j,1) = exp(sh%lfac(m) + rm*sh%logs(j))   ! underflow -> 0 is harmless
    end do
    if (nl >= 2) then
      a = sqrt(real(2*m+3,REAL64))
      P(:,2) = a * sh%x(:) * P(:,1)
    end if
    do k = 3, nl
      l = m + k - 1
      rl = real(l,REAL64)
      a = sqrt((4.0d0*rl*rl - 1.0d0) / (rl*rl - rm*rm))
      b = sqrt(((rl-1.0d0)**2 - rm*rm) / (4.0d0*(rl-1.0d0)**2 - 1.0d0))
      do j = 1, sh%nh
        P(j,k) = a * (sh%x(j)*P(j,k-1) - b*P(j,k-2))
      end do
    end do
  end subroutine plm_block

  ! ======================================================================
  ! sh_draw: two independent realisations f1, f2 on the (nlon, nlat) grid
  ! (lat index 1 = southernmost ring).  Every rank returns the full fields.

  subroutine sh_draw(sh, rseed, f1, f2)
    type(sphere_harmonics_random_fields), intent(inout) :: sh
    integer, intent(inout) :: rseed(NRANDSEED)
    real, intent(out) :: f1(sh%nphi_sphere, sh%nz_sphere), f2(sh%nphi_sphere, sh%nz_sphere)

    complex(REAL64), allocatable :: G(:,:,:), Gl(:)
    real(REAL64), allocatable :: Pw(:), Be(:,:), Bo(:,:), E(:,:), O(:,:), gz(:,:), rout(:)
    real, allocatable :: u(:), loc1(:), loc2(:)
    integer, allocatable :: rc(:), rd(:)
    integer :: i, m, nl, ne, no, k, l, j, c, ierr, nlat, nlon, nu
    integer(INT64) :: n0
    real(REAL64) :: fac, r1, r2

    nlat = sh%nz_sphere; nlon = sh%nphi_sphere

    ! --- global uniform stream: 4 normals per (l,m), Box-Muller pairs ----
    nu = int(4_INT64*sh%ncoef)
    allocate(u(nu))
    call nr_ran2_2d(1, nu, rseed, u)

    allocate(G(sh%nm, 2, nlat)); G = (0.0d0, 0.0d0)
    allocate(Be((sh%L+2)/2, 4), Bo((sh%L+2)/2, 4), E(sh%nh,4), O(sh%nh,4))
    allocate(gz(sh%L+1, 4))
    if (.not. sh%use_table) allocate(Pw(int(sh%nh,INT64)*int(sh%L+1,INT64)))

    do i = 1, sh%nmloc
      m  = sh%mloc(i)
      nl = sh%L - m + 1
      ne = (nl+1)/2; no = nl/2
      fac = merge(1.0d0, sqrt(0.5d0), m == 0)
      do k = 1, nl
        l  = m + k - 1
        n0 = 4_INT64*(sh%coff(m) + k - 1)
        do c = 1, 3, 2
          r1 = sqrt(-2.0d0*log(max(real(u(n0+c),REAL64), 1.0d-300)))
          r2 = 2.0d0*PI*real(u(n0+c+1),REAL64)
          gz(k,c)   = r1*cos(r2)
          gz(k,c+1) = r1*sin(r2)
        end do
        gz(k,:) = gz(k,:) * fac * sh%sqrtC(l) * sh%scale
      end do
      if (m == 0) then
        gz(1:nl,2) = 0.0d0; gz(1:nl,4) = 0.0d0
      end if
      Be(1:ne,:) = gz(1:nl:2,:)
      if (no > 0) Bo(1:no,:) = gz(2:nl:2,:)

      if (sh%use_table) then
        call dgemm('N','N', sh%nh, 4, ne, 1.0d0, sh%ptab(sh%off(i)), 2*sh%nh, &
                   Be, size(Be,1), 0.0d0, E, sh%nh)
        if (no > 0) then
          call dgemm('N','N', sh%nh, 4, no, 1.0d0, sh%ptab(sh%off(i)+sh%nh), 2*sh%nh, &
                     Bo, size(Bo,1), 0.0d0, O, sh%nh)
        else
          O = 0.0d0
        end if
      else
        call plm_block(sh, m, nl, Pw(1:int(sh%nh,INT64)*nl))
        call dgemm('N','N', sh%nh, 4, ne, 1.0d0, Pw(1), 2*sh%nh, &
                   Be, size(Be,1), 0.0d0, E, sh%nh)
        if (no > 0) then
          call dgemm('N','N', sh%nh, 4, no, 1.0d0, Pw(1+sh%nh), 2*sh%nh, &
                     Bo, size(Bo,1), 0.0d0, O, sh%nh)
        else
          O = 0.0d0
        end if
      end if

      ! north ring j (theta_j) -> lat index nlat+1-j ; south mirror -> j
      do j = 1, sh%nh
        G(m+1,1,nlat+1-j) = cmplx(E(j,1)+O(j,1), E(j,2)+O(j,2), REAL64)
        G(m+1,2,nlat+1-j) = cmplx(E(j,3)+O(j,3), E(j,4)+O(j,4), REAL64)
        G(m+1,1,j)        = cmplx(E(j,1)-O(j,1), E(j,2)-O(j,2), REAL64)
        G(m+1,2,j)        = cmplx(E(j,3)-O(j,3), E(j,4)-O(j,4), REAL64)
      end do
    end do
    deallocate(u)

    ! --- transpose m-distributed -> ring-distributed --------------------
    allocate(Gl(max(1, 2*sh%nm*sh%nloc)))
    if (sh%comm /= MPI_COMM_NULL) then
      allocate(rc(sh%npes))
      rc = 2*sh%nm*sh%rcount
      call MPI_Reduce_scatter(G, Gl, rc, MPI_DOUBLE_COMPLEX, MPI_SUM, sh%comm, ierr)
    else
      Gl = reshape(G, [2*sh%nm*nlat])
    end if
    deallocate(G)

    ! --- inverse real FFT of local rings (both fields) ------------------
    allocate(rout(max(1, nlon*2*sh%nloc)))
#ifdef MKL_AVAILABLE
    if (sh%nloc > 0) ierr = DftiComputeBackward(sh%desc, Gl, rout)
#else
    call irfft_rows(Gl, 2*sh%nloc, sh%nm, nlon, rout)
#endif

    if (sh%comm /= MPI_COMM_NULL) then
      allocate(loc1(max(1,nlon*sh%nloc)), loc2(max(1,nlon*sh%nloc)))
      do j = 1, sh%nloc
        loc1((j-1)*nlon+1:j*nlon) = real(rout((2*j-2)*nlon+1:(2*j-1)*nlon))
        loc2((j-1)*nlon+1:j*nlon) = real(rout((2*j-1)*nlon+1:(2*j)*nlon))
      end do
      allocate(rd(sh%npes))
      rc = nlon*sh%rcount; rd = nlon*sh%rdispl
      call MPI_Allgatherv(loc1, nlon*sh%nloc, MPI_REAL, f1, rc, rd, MPI_REAL, sh%comm, ierr)
      call MPI_Allgatherv(loc2, nlon*sh%nloc, MPI_REAL, f2, rc, rd, MPI_REAL, sh%comm, ierr)
    else
      do j = 1, nlat
        f1(:,j) = real(rout((2*j-2)*nlon+1:(2*j-1)*nlon))
        f2(:,j) = real(rout((2*j-1)*nlon+1:(2*j)*nlon))
      end do
    end if
  end subroutine sh_draw

#ifndef MKL_AVAILABLE
  ! Pure-Fortran O(nph^2) inverse real DFT without 1/N (matches MKL backward).
  subroutine irfft_rows(inp, nrows, nm, nph, out)
    integer, intent(in) :: nrows, nm, nph
    complex(REAL64), intent(in) :: inp(nm, *)
    real(REAL64), intent(out) :: out(nph, *)
    integer :: row, j, k
    real(REAL64) :: ang, s, w
    w = 2.0d0*PI/real(nph,REAL64)
    do row = 1, nrows
      do j = 0, nph-1
        s = real(inp(1,row))
        do k = 1, nm-2
          ang = w*real(j,REAL64)*real(k,REAL64)
          s = s + 2.0d0*(real(inp(k+1,row))*cos(ang) - aimag(inp(k+1,row))*sin(ang))
        end do
        ang = w*real(j,REAL64)*real(nm-1,REAL64)
        s = s + real(inp(nm,row))*cos(ang)
        out(j+1,row) = s
      end do
    end do
  end subroutine irfft_rows
#endif

  ! ======================================================================
  ! Copy, or bilinearly remap (periodic in lon, cell-centred lat) the native
  ! (nphi, nz) field onto an arbitrary (N_lon, N_lat) cell-centred grid.

  subroutine remap_latlon(src, dst)
    real, intent(in)  :: src(:,:)
    real, intent(out) :: dst(:,:)
    integer :: ns, nzs, nd, nzd, i, j, p1, p2, k1, k2
    real(REAL64) :: pos, wp, wz, v1, v2

    ns = size(src,1); nzs = size(src,2)
    nd = size(dst,1); nzd = size(dst,2)
    if (ns == nd .and. nzs == nzd) then
      dst = src
      return
    end if
    do j = 1, nzd
      pos = 0.5d0 + (real(j,REAL64) - 0.5d0) * real(nzs,REAL64) / real(nzd,REAL64)
      k1 = max(1, min(nzs, floor(pos)))
      k2 = max(1, min(nzs, k1 + 1))
      wz = max(0.0d0, min(1.0d0, pos - real(k1,REAL64)))
      do i = 1, nd
        pos = 0.5d0 + (real(i,REAL64) - 0.5d0) * real(ns,REAL64) / real(nd,REAL64)
        p1 = floor(pos)
        wp = pos - real(p1,REAL64)
        p1 = modulo(p1 - 1, ns) + 1
        p2 = modulo(p1, ns) + 1
        v1 = (1.0d0-wz)*src(p1,k1) + wz*src(p1,k2)
        v2 = (1.0d0-wz)*src(p2,k1) + wz*src(p2,k2)
        dst(i,j) = real((1.0d0-wp)*v1 + wp*v2)
      end do
    end do
  end subroutine remap_latlon

  subroutine normalize_unit_variance(field)
    real, intent(inout) :: field(:,:)
    real(REAL64) :: field_mean, field_norm
    ! Accumulate in double precision: single-precision sums over ~4e6
    ! points bias the variance by O(1e-3).
    field_mean = sum(real(field,REAL64)) / real(size(field),REAL64)
    field = real(real(field,REAL64) - field_mean)
    field_norm = sqrt(sum(real(field,REAL64)**2) / real(size(field),REAL64))
    if (field_norm > 0.0_REAL64) field = real(real(field,REAL64) / field_norm)
  end subroutine normalize_unit_variance

  ! ======================================================================
  ! xcorr_deg_to_c2: calibrate c2 from xcorr.
  ! xcorr is the Gaussian sigma: C(xcorr)/C(0) = exp(-0.5).

  subroutine xcorr_deg_to_c2(xcorr_deg, nz, c2_out)
    real(REAL64), intent(in)  :: xcorr_deg
    integer,        intent(in)  :: nz
    real(REAL64), intent(out) :: c2_out
    real(REAL64) :: lo, hi, mid, r, z, target
    integer :: it, lmax

    z = cos(xcorr_deg*MAPL_DEGREES_TO_RADIANS_R8)
    lmax = 3*nz
    target = exp(-0.5d0)
    lo = 1.0d-12; hi = 1.0d2
    if (ratio(lo) > target) then
      c2_out = lo; return
    end if
    if (ratio(hi) < target) then
      c2_out = hi; return
    end if
    mid = lo
    do it = 1, 80
      mid = 0.5d0*(lo+hi)
      r = ratio(mid)
      if (abs(r - target) < 1.0d-12) exit
      if (r > target) then
        hi = mid
      else
        lo = mid
      end if
    end do
    c2_out = mid
  contains
    real(REAL64) function ratio(c2)
      real(REAL64), intent(in) :: c2
      real(REAL64) :: num, den, C, P0, P1, P2
      integer :: l
      num = 0; den = 0; P0 = 1; P1 = z
      do l = 0, lmax
        C = 1.0d0/(1.0d0 + c2*(real(l,REAL64)*real(l+1,REAL64))**2)
        den = den + real(2*l+1,REAL64)*C
        if (l == 0) then
          num = num + C
        else if (l == 1) then
          num = num + 3.0d0*C*z
        else
          P2 = (real(2*l-1,REAL64)*z*P1 - real(l-1,REAL64)*P0)/real(l,REAL64)
          num = num + real(2*l+1,REAL64)*C*P2
          P0 = P1; P1 = P2
        end if
      end do
      ratio = num/den
    end function ratio
  end subroutine xcorr_deg_to_c2

end module sphere_harmonics_random_fields_mod
