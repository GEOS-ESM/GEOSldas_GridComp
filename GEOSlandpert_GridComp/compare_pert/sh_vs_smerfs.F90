! Comparison: SMERFS sphere filter (sphere_random_fields) vs direct
! spherical-harmonic synthesis (sphere_harmonics_random_fields).
!
!   mpirun -np N sh_vs_smerfs.x [xcorr] [nlon] [nlat] [ndraw]
!
! Both generators are driven through the same abstract interface
! (constructor + generate_2d_Random_field with an NR rseed).  Reported:
!   * build time, first/min/mean draw time
!   * field mean / std
!   * empirical correlation vs. analytic Legendre-series correlation for
!     lags of 0.5,1,2,3 x xcorr along longitude (equator, 60N) and along
!     latitude (equator), averaged over all draws and both fields.
program sh_vs_smerfs
  use mpi
  use nr_ran2_gasdev, only: NRANDSEED, init_randseed
  use abstract_random_fieldsMod, only: abstract_random_fields
  use sphere_random_fields_mod, only: sphere_random_fields
  use sphere_harmonics_random_fields_mod, only: sphere_harmonics_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  implicit none

  integer, parameter :: NLAG = 4, NCASE = 3, MAXDRAW = 50
  real(8), parameter :: PI = 3.14159265358979323846d0
  real(8), parameter :: lagmul(NLAG) = [0.5d0, 1.0d0, 2.0d0, 3.0d0]
  character(len=*), parameter :: casename(NCASE) = &
       [character(len=10) :: 'lon @ eq', 'lon @ 60N', 'lat @ eq']

  integer :: ierr, rank, npes, nlon, nlat, ndraw
  real :: xcorr
  character(len=32) :: arg
  type(pert_param_type) :: pp
  type(grid_def_type) :: grid
  real, allocatable :: f1(:,:), f2(:,:)
  real(8) :: corr_sm(NLAG,NCASE), corr_sh(NLAG,NCASE), corr_sht(NLAG,NCASE)
  real(8) :: c2_sm, c2_sh

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
  xcorr = 1.0; nlon = 2880; nlat = 1440; ndraw = 3
  call get_command_argument(1, arg); if (len_trim(arg) > 0) read(arg,*) xcorr
  call get_command_argument(2, arg); if (len_trim(arg) > 0) read(arg,*) nlon
  call get_command_argument(3, arg); if (len_trim(arg) > 0) read(arg,*) nlat
  call get_command_argument(4, arg); if (len_trim(arg) > 0) read(arg,*) ndraw
  ndraw = max(1, min(ndraw, MAXDRAW))

  pp%xcorr = xcorr; pp%ycorr = xcorr
  grid%N_lon = nlon; grid%N_lat = nlat
  grid%dlon = 360.0/nlon; grid%dlat = 180.0/nlat
  allocate(f1(nlon,nlat), f2(nlon,nlat))

  if (rank == 0) write(*,'(a,i0,a,i0,a,i0,a,f6.3,a,i0)') 'ranks=', npes, ' grid=', nlon, 'x', nlat, &
       ' xcorr=', xcorr, ' ndraw=', ndraw

  call run('SMERFS        ', 1, corr_sm, c2_sm)
  call run('SH table      ', 2, corr_sh,  c2_sh)
  call run('SH on-the-fly ', 3, corr_sht, c2_sh)

  if (rank == 0) then
    write(*,'(/,a,es12.5,a,es12.5)') 'c2: SMERFS=', c2_sm, '  SH=', c2_sh
    write(*,'(a)') 'Correlation vs lag (analytic from c2, empirical averaged over draws):'
    write(*,'(a10,a8,4a11)') 'case', 'lag/xc', 'analytic', 'SMERFS', 'SH table', 'SH otf'
    block
      integer :: ic, il
      do ic = 1, NCASE
        do il = 1, NLAG
          write(*,'(a10,f8.2,4f11.4)') casename(ic), lagmul(il), &
               analytic_corr(ic, il, c2_sh), corr_sm(il,ic), corr_sh(il,ic), corr_sht(il,ic)
        end do
      end do
    end block
    write(*,'(a,f8.5)') 'max |SH table - SH otf| corr diff (should be ~0): ', &
         maxval(abs(corr_sh - corr_sht))
  end if

  call MPI_Finalize(ierr)

contains

  ! kind: 1 = SMERFS, 2 = SH with Legendre table, 3 = SH on-the-fly
  subroutine run(tag, kind, corr, c2)
    character(*), intent(in) :: tag
    integer, intent(in) :: kind
    real(8), intent(out) :: corr(NLAG,NCASE), c2
    class(abstract_random_fields), allocatable :: rf
    integer :: rseed(NRANDSEED), k
    real(8) :: t0, tb, td(MAXDRAW), mean, sd, msum, ssum
    integer(8) :: tbytes

    rseed = 0; rseed(1) = -777; call init_randseed(rseed)
    tbytes = 0

    call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = MPI_Wtime()
    select case (kind)
    case (1)
      allocate(rf, source=sphere_random_fields(pp, grid, comm=MPI_COMM_WORLD))
    case (2)
      allocate(rf, source=sphere_harmonics_random_fields(pp, grid, comm=MPI_COMM_WORLD, use_table=.true.))
    case default
      allocate(rf, source=sphere_harmonics_random_fields(pp, grid, comm=MPI_COMM_WORLD, use_table=.false.))
    end select
    call MPI_Barrier(MPI_COMM_WORLD, ierr); tb = MPI_Wtime() - t0

    select type (rf)
    type is (sphere_random_fields)
      c2 = rf%c2
    type is (sphere_harmonics_random_fields)
      c2 = rf%c2; tbytes = rf%table_bytes
    end select

    corr = 0; msum = 0; ssum = 0
    do k = 1, ndraw
      call MPI_Barrier(MPI_COMM_WORLD, ierr); t0 = MPI_Wtime()
      call rf%generate_2d_Random_field(rseed, f1, f2, xcorr, xcorr, grid%dlon, grid%dlat)
      call MPI_Barrier(MPI_COMM_WORLD, ierr); td(k) = MPI_Wtime() - t0
      if (rank == 0) then
        call accumulate(f1, corr, msum, ssum)
        call accumulate(f2, corr, msum, ssum)
      end if
    end do
    corr = corr / (2*ndraw)
    mean = msum / (2*ndraw); sd = ssum / (2*ndraw)

    if (rank == 0) then
      write(*,'(a,a,f10.4,a,f9.4,a,f9.4,a,f9.4,a,f8.4,a,f8.4)') tag, ' build=', tb, &
        ' draw1=', td(1), ' drawmin=', minval(td(1:ndraw)), ' drawmean=', sum(td(1:ndraw))/ndraw, &
        ' | mean=', mean, ' std=', sd
      if (kind == 2) write(*,'(a,f8.3,a)') '  SH Legendre table total = ', &
           real(tbytes,8)/1024.0d0**3, ' GiB'
    end if
    call rf%finalize()
  end subroutine run

  subroutine accumulate(f, corr, msum, ssum)
    real, intent(in) :: f(:,:)
    real(8), intent(inout) :: corr(NLAG,NCASE), msum, ssum
    real(8) :: m
    integer :: il, ic
    m = sum(real(f,8))/size(f)
    msum = msum + m
    ssum = ssum + sqrt(sum((real(f,8)-m)**2)/(size(f)-1))
    do ic = 1, NCASE
      do il = 1, NLAG
        corr(il,ic) = corr(il,ic) + empirical_corr(f, ic, il)
      end do
    end do
  end subroutine accumulate

  ! Band of rings around the target latitude (+-5 deg), lag in grid cells.
  real(8) function empirical_corr(f, ic, il)
    real, intent(in) :: f(:,:)
    integer, intent(in) :: ic, il
    integer :: i, j, j0, jb, lag, n
    real(8) :: s, a, b
    if (ic == 3) then
      lag = max(1, nint(lagmul(il)*xcorr / grid%dlat))
      j0 = nlat/2
    else
      lag = max(1, nint(lagmul(il)*xcorr / grid%dlon))
      j0 = merge(nlat/2, nint((60.0+90.0)/grid%dlat), ic == 1)
    end if
    jb = max(1, nint(5.0/grid%dlat))
    s = 0; a = 0; b = 0; n = 0
    do j = max(1, j0-jb), min(nlat - merge(lag, 0, ic == 3), j0+jb)
      do i = 1, nlon
        if (ic == 3) then
          s = s + real(f(i,j),8)*f(i,j+lag)
          b = b + real(f(i,j+lag),8)**2
        else
          s = s + real(f(i,j),8)*f(mod(i-1+lag,nlon)+1,j)
          b = b + real(f(mod(i-1+lag,nlon)+1,j),8)**2
        end if
        a = a + real(f(i,j),8)**2
      end do
    end do
    empirical_corr = s / sqrt(a*b)
  end function empirical_corr

  ! Analytic C(theta)/C(0) for C_l = 1/(1 + c2 (l(l+1))^2), at the
  ! great-circle angle of the grid lag used by empirical_corr.
  real(8) function analytic_corr(ic, il, c2)
    integer, intent(in) :: ic, il
    real(8), intent(in) :: c2
    real(8) :: z, lat, dl, num, den, C, P0, P1, P2
    integer :: l, lag
    if (ic == 3) then
      lag = max(1, nint(lagmul(il)*xcorr / grid%dlat))
      z = cos(lag*grid%dlat*PI/180d0)
    else
      lag = max(1, nint(lagmul(il)*xcorr / grid%dlon))
      lat = merge(0.0d0, 60.0d0, ic == 1)*PI/180d0
      dl  = lag*grid%dlon*PI/180d0
      z = sin(lat)**2 + cos(lat)**2*cos(dl)
    end if
    num = 0; den = 0; P0 = 1; P1 = z
    do l = 0, 3*nlat
      C = 1.0d0/(1.0d0 + c2*(real(l,8)*real(l+1,8))**2)
      den = den + real(2*l+1,8)*C
      if (l == 0) then
        num = num + C
      else if (l == 1) then
        num = num + 3.0d0*C*z
      else
        P2 = (real(2*l-1,8)*z*P1 - real(l-1,8)*P0)/real(l,8)
        num = num + real(2*l+1,8)*C*P2
        P0 = P1; P1 = P2
      end if
    end do
    analytic_corr = num/den
  end function analytic_corr

end program sh_vs_smerfs
