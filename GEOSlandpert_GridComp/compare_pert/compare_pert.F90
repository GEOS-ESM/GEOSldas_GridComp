! compare_pert: performance and output comparison of
!   rectangle_random_fields (flat lat-lon FFT) and
!   sphere_random_fields    (spherical-harmonic synthesis)
! on a global regular lat-lon grid.
!
! Usage:
!   mpirun -np N compare_pert.x [xcorr_deg] [ndraw] [nlon] [nlat]
!
! For each generator:
!   - construction (init) time, max over ranks
!   - ndraw calls to generate_2d_Random_field (2 fields per call);
!     mean and min wall time per call
! Rank 0 writes the first field of the first draw to latlon_pert.bin /
! sphere_pert.bin (real(4) stream, lon fastest) with a small header file
! compare_pert_grid.txt, and appends one line per generator to
! compare_pert_timing.csv:
!   generator,npes,nlon,nlat,xcorr,ndraw,init_s,draw_mean_s,draw_min_s

program compare_pert
  use mpi
  use nr_ran2_gasdev, only: NRANDSEED, init_randseed
  use abstract_random_fieldsMod, only: abstract_random_fields
  use rectangle_random_fieldsMod, only: rectangle_random_fields
  use sphere_random_fields_mod, only: sphere_harmonics_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  implicit none

  integer, parameter :: RSEED0 = -777
  integer :: mpierr, rank, npes, ndraw, nlon, nlat, unit
  integer :: seed_template(NRANDSEED)
  real :: xcorr, dlon, dlat
  character(len=32) :: arg

  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)

  xcorr = 2.0; ndraw = 10; nlon = 720; nlat = 360
  call get_command_argument(1, arg); if (len_trim(arg) > 0) read(arg,*) xcorr
  call get_command_argument(2, arg); if (len_trim(arg) > 0) read(arg,*) ndraw
  call get_command_argument(3, arg); if (len_trim(arg) > 0) read(arg,*) nlon
  call get_command_argument(4, arg); if (len_trim(arg) > 0) read(arg,*) nlat
  if (xcorr <= 0.0 .or. ndraw < 1 .or. nlon < 2 .or. nlat < 2) then
    if (rank == 0) write(*,'(a)') 'Usage: compare_pert.x [xcorr_deg>0] [ndraw>=1] [nlon] [nlat]'
    call MPI_Abort(MPI_COMM_WORLD, 1, mpierr)
  end if
  dlon = 360.0 / real(nlon); dlat = 180.0 / real(nlat)

  if (rank == 0) then
    write(*,'(a,i0,a,i0,a,i0,a,f7.3,a,i0)') '=== compare_pert: ranks=', npes, &
         ' grid=', nlon, 'x', nlat, ' xcorr=', xcorr, ' ndraw=', ndraw
    open(newunit=unit, file='compare_pert_grid.txt', status='replace')
    write(unit,*) nlon, nlat, xcorr
    close(unit)
  end if

  ! One RNG state, an identical copy given to each generator.
  seed_template = 0
  seed_template(1) = RSEED0
  call init_randseed(seed_template)

  call run(1, 'rectangle', 'latlon_pert.bin')
  call run(2, 'sphere',    'sphere_pert.bin')

  call MPI_Finalize(mpierr)

contains

  subroutine run(kind, tag, fname)
    integer, intent(in) :: kind
    character(*), intent(in) :: tag, fname
    class(abstract_random_fields), allocatable :: gen
    type(pert_param_type) :: pp
    type(grid_def_type) :: grid
    integer :: rseed(NRANDSEED), i, u
    real, allocatable :: f1(:,:), f2(:,:)
    double precision :: t0, t, t_init, t_sum, t_min
    logical :: exists

    pp%xcorr = xcorr; pp%ycorr = xcorr; pp%coarsen = .false.
    grid%N_lon = nlon; grid%N_lat = nlat
    grid%dlon = dlon;  grid%dlat = dlat
    rseed = seed_template
    allocate(f1(nlon,nlat), f2(nlon,nlat))

    call MPI_Barrier(MPI_COMM_WORLD, mpierr); t0 = MPI_Wtime()
    select case (kind)
    case (1)
#ifdef MKL_AVAILABLE
      allocate(gen, source=rectangle_random_fields(pp, grid, comm=MPI_COMM_WORLD))
#else
      allocate(gen, source=rectangle_random_fields(pp, grid))
#endif
    case (2)
      allocate(gen, source=sphere_harmonics_random_fields(pp, grid, comm=MPI_COMM_WORLD))
    end select
    call MPI_Barrier(MPI_COMM_WORLD, mpierr); t_init = MPI_Wtime() - t0

    t_sum = 0d0; t_min = huge(1d0)
    do i = 1, ndraw
      call MPI_Barrier(MPI_COMM_WORLD, mpierr); t0 = MPI_Wtime()
      call gen%generate_2d_Random_field(rseed, f1, f2, xcorr, xcorr, dlon, dlat)
      call MPI_Barrier(MPI_COMM_WORLD, mpierr); t = MPI_Wtime() - t0
      t_sum = t_sum + t; t_min = min(t_min, t)
      if (i == 1 .and. rank == 0) then
        open(newunit=u, file=fname, status='replace', access='stream', form='unformatted')
        write(u) f1; close(u)
        call print_stats(tag, f1)
      end if
    end do

    if (rank == 0) then
      write(*,'(a,f10.4,a,f10.4,a,f10.4)') '  init[s]=', t_init, &
           '  draw_mean[s]=', t_sum/ndraw, '  draw_min[s]=', t_min
      inquire(file='compare_pert_timing.csv', exist=exists)
      if (exists) then
        open(newunit=u, file='compare_pert_timing.csv', status='old', position='append')
      else
        open(newunit=u, file='compare_pert_timing.csv', status='new')
        write(u,'(a)') 'generator,npes,nlon,nlat,xcorr,ndraw,init_s,draw_mean_s,draw_min_s'
      end if
      write(u,'(a,",",i0,",",i0,",",i0,",",f0.4,",",i0,3(",",es12.5))') &
           tag, npes, nlon, nlat, xcorr, ndraw, t_init, t_sum/ndraw, t_min
      close(u)
    end if
    call gen%finalize()
  end subroutine run

  subroutine print_stats(tag, field)
    character(*), intent(in) :: tag
    real, intent(in) :: field(:,:)
    real :: mean, std
    mean = sum(field) / real(size(field))
    std = sqrt(sum((field-mean)**2) / real(size(field)-1))
    write(*,'(a,a10,a,f9.5,a,f9.5,a,f9.5,a,f9.5)') '  [',tag,'] mean=',mean, &
      ' std=',std,' min=',minval(field),' max=',maxval(field)
  end subroutine print_stats

end program compare_pert
