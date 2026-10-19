! Timing comparison: flat lat-lon FFT generator (rectangle_random_fields)
! vs direct spherical-harmonic synthesis (sphere_harmonics_random_fields).
!
!   mpirun -np N sh_vs_rectangle.x [xcorr] [nlon] [nlat] [ndraw]
!
! Cases (all driven through the abstract_random_fields interface, same rseed):
!   RECT full      - rectangle_random_fields, coarsen = .false.
!   RECT coarsen   - rectangle_random_fields, coarsen = .true. (output on the
!                    coarse N_x x N_y grid, as used in land_pert)
!   SH             - sphere_harmonics_random_fields (default: no table)
! Reported: build time, first/min/mean draw time (2 fields per draw),
! field mean/std and correlation at lag = xcorr along longitude at the
! equator and at 60N, and along latitude at the equator, plus the
! analytic sphere value for reference.  The rectangle generator works in
! degrees on a flat grid, so its longitudinal correlation length in km
! shrinks with cos(lat); the SH field is isotropic on the sphere.
program sh_vs_rectangle
  use mpi
  use, intrinsic :: iso_fortran_env, only: REAL64
  use nr_ran2_gasdev, only: NRANDSEED, init_randseed
  use abstract_random_fieldsMod, only: abstract_random_fields
  use rectangle_random_fieldsMod, only: rectangle_random_fields
  use sphere_random_fields_mod, only: sphere_harmonics_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  implicit none

  integer, parameter :: MAXDRAW = 50
  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64

  integer :: ierr, rank, npes, nlon, nlat, ndraw
  real :: xcorr
  real(REAL64) :: c2_sh
  character(len=32) :: arg
  type(grid_def_type) :: grid

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
  xcorr = 1.0; nlon = 2880; nlat = 1440; ndraw = 10
  call get_command_argument(1, arg); if (len_trim(arg) > 0) read(arg,*) xcorr
  call get_command_argument(2, arg); if (len_trim(arg) > 0) read(arg,*) nlon
  call get_command_argument(3, arg); if (len_trim(arg) > 0) read(arg,*) nlat
  call get_command_argument(4, arg); if (len_trim(arg) > 0) read(arg,*) ndraw
  ndraw = max(1, min(ndraw, MAXDRAW))
  grid%N_lon = nlon; grid%N_lat = nlat
  grid%dlon = 360.0/nlon; grid%dlat = 180.0/nlat

  if (rank == 0) then
    write(*,'(a,i0,a,i0,a,i0,a,f6.3,a,i0)') 'ranks=', npes, ' grid=', nlon, 'x', nlat, &
         ' xcorr=', xcorr, ' ndraw=', ndraw
#ifndef MKL_AVAILABLE
    write(*,'(a)') 'NOTE: built without MKL; rectangle runs serially on every rank'
#endif
    write(*,'(a14,a10,4a10,a9,a8,3a9)') 'case', 'grid', 'build', 'draw1', 'drawmean', &
         'drawmin', 'std', 'lon@eq', 'lon@60N', 'lat@eq'
  end if

  c2_sh = 0
  call run('SH (no table)', 3)
  call run('RECT full', 1)
  call run('RECT coarsen', 2)

  if (rank == 0) write(*,'(a,3f9.4,a)') 'analytic sphere corr at lag xcorr:       ', &
       analytic_corr(1), analytic_corr(2), analytic_corr(3), '  (SH spectrum)'

  call MPI_Finalize(ierr)

contains

  ! kind: 1 = rectangle full, 2 = rectangle coarsened, 3 = SH
  subroutine run(tag, kind)
    character(*), intent(in) :: tag
    integer, intent(in) :: kind
    class(abstract_random_fields), allocatable :: rf
    type(pert_param_type) :: pp
    real, allocatable :: f1(:,:), f2(:,:)
    integer :: rseed(NRANDSEED), k, nx, ny
    real :: dx, dy
    real(REAL64) :: t0, tb, td(MAXDRAW), sd, c(3)
    character(len=16) :: gstr

    pp%xcorr = xcorr; pp%ycorr = xcorr
    pp%coarsen = (kind == 2)
    rseed = 0; rseed(1) = -777; call init_randseed(rseed)

    call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = MPI_Wtime()
    select case (kind)
    case (1, 2)
#ifdef MKL_AVAILABLE
      allocate(rf, source=rectangle_random_fields(pp, grid, comm=MPI_COMM_WORLD))
#else
      allocate(rf, source=rectangle_random_fields(pp, grid))
#endif
    case default
      allocate(rf, source=sphere_harmonics_random_fields(pp, grid, comm=MPI_COMM_WORLD))
    end select
    call MPI_Barrier(MPI_COMM_WORLD, ierr); tb = MPI_Wtime() - t0

    nx = nlon; ny = nlat; dx = grid%dlon; dy = grid%dlat
    select type (rf)
    type is (rectangle_random_fields)
      nx = rf%N_x; ny = rf%N_y; dx = rf%rdlon; dy = rf%rdlat
    type is (sphere_harmonics_random_fields)
      c2_sh = rf%c2
    end select
    allocate(f1(nx,ny), f2(nx,ny))

    sd = 0; c = 0
    do k = 1, ndraw
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = MPI_Wtime()
      call rf%generate_2d_Random_field(rseed, f1, f2, xcorr, xcorr, dx, dy)
      call MPI_Barrier(MPI_COMM_WORLD, ierr); td(k) = MPI_Wtime() - t0
      if (rank == 0 .and. k == 1) call dump(kind, f1)
      if (rank == 0) then
        sd = sd + std(f1) + std(f2)
        c(1) = c(1) + lagcorr(f1, dx, dy, 0.0, .false.) + lagcorr(f2, dx, dy, 0.0, .false.)
        c(2) = c(2) + lagcorr(f1, dx, dy, 60.0, .false.) + lagcorr(f2, dx, dy, 60.0, .false.)
        c(3) = c(3) + lagcorr(f1, dx, dy, 0.0, .true.) + lagcorr(f2, dx, dy, 0.0, .true.)
      end if
    end do
    sd = sd/(2*ndraw); c = c/(2*ndraw)

    if (rank == 0) then
      write(gstr,'(i0,"x",i0)') nx, ny
      write(*,'(a14,a10,4f10.4,f9.4,f8.4,3f9.4)') tag, trim(gstr), tb, td(1), &
           sum(td(1:ndraw))/ndraw, minval(td(1:ndraw)), sd, c
    end if
    call rf%finalize()
  end subroutine run

  ! First field of the first draw -> stream binary: int32 nx, ny, then
  ! real(4) f(nx,ny) (lon fastest, lat from south).  Read by plot_sh_vs_rectangle.py.
  subroutine dump(kind, f)
    integer, intent(in) :: kind
    real, intent(in) :: f(:,:)
    character(len=*), parameter :: fname(3) = &
         [character(len=24) :: 'rect_pert.bin', 'rect_coarsen_pert.bin', 'sh_pert.bin']
    integer :: u
    open(newunit=u, file=trim(fname(kind)), status='replace', access='stream', form='unformatted')
    write(u) size(f,1), size(f,2), f
    close(u)
  end subroutine dump

  real(REAL64) function std(f)
    real, intent(in) :: f(:,:)
    real(REAL64) :: m
    m = sum(real(f,REAL64))/size(f)
    std = sqrt(sum((real(f,REAL64)-m)**2)/(size(f)-1))
  end function std

  ! Correlation at lag = xcorr (in grid cells of spacing dx/dy) in a +-5 deg
  ! band around latitude lat0; along longitude (periodic) or latitude.
  real(REAL64) function lagcorr(f, dx, dy, lat0, along_lat)
    real, intent(in) :: f(:,:), dx, dy, lat0
    logical, intent(in) :: along_lat
    integer :: i, j, nx, ny, j0, jb, lag, jhi
    real(REAL64) :: s, a, b, v0, v1
    nx = size(f,1); ny = size(f,2)
    if (along_lat) then
      lag = max(1, nint(xcorr/dy))
    else
      lag = max(1, nint(xcorr/dx))
    end if
    j0 = max(1, min(ny, nint((lat0 + 90.0)/dy + 0.5)))
    jb = max(1, nint(5.0/dy))
    jhi = min(ny, j0+jb)
    if (along_lat) jhi = min(ny-lag, j0+jb)
    s = 0; a = 0; b = 0
    do j = max(1, j0-jb), jhi
      do i = 1, nx
        v0 = f(i,j)
        if (along_lat) then
          v1 = f(i,j+lag)
        else
          v1 = f(mod(i-1+lag,nx)+1,j)
        end if
        s = s + v0*v1; a = a + v0*v0; b = b + v1*v1
      end do
    end do
    lagcorr = s/sqrt(a*b)
  end function lagcorr

  ! Analytic C(theta)/C(0) of the SH spectrum at lag xcorr (grid-rounded).
  real(REAL64) function analytic_corr(icase)
    integer, intent(in) :: icase
    real(REAL64) :: z, lat, dl, num, den, C, P0, P1, P2
    integer :: l, lag
    if (icase == 3) then
      lag = max(1, nint(xcorr/grid%dlat))
      z = cos(lag*grid%dlat*PI/180)
    else
      lag = max(1, nint(xcorr/grid%dlon))
      lat = merge(0.0_REAL64, 60.0_REAL64, icase == 1)*PI/180
      dl  = lag*grid%dlon*PI/180
      z = sin(lat)**2 + cos(lat)**2*cos(dl)
    end if
    num = 0; den = 0; P0 = 1; P1 = z
    do l = 0, 3*nlat
      C = 1/(1 + c2_sh*(real(l,REAL64)*real(l+1,REAL64))**2)
      den = den + (2*l+1)*C
      if (l == 0) then
        num = num + C
      else if (l == 1) then
        num = num + 3*C*z
      else
        P2 = ((2*l-1)*z*P1 - (l-1)*P0)/l
        num = num + (2*l+1)*C*P2
        P0 = P1; P1 = P2
      end if
    end do
    analytic_corr = num/den
  end function analytic_corr

end program sh_vs_rectangle
