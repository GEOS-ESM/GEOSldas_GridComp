program compare_pert
  use mpi
  use, intrinsic :: iso_c_binding, only: c_double
  use nr_ran2_gasdev, only: NRANDSEED, init_randseed
  use rectangle_random_fieldsMod, only: rectangle_random_fields
  use sphere_random_fields_mod, only: sphere_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  implicit none

  integer, parameter :: N_LON = 720, N_LAT = 360
  real, parameter :: DLON = 360.0 / real(N_LON), DLAT = 180.0 / real(N_LAT)
  real, parameter :: XCORR = 10.0, YCORR = 10.0
  integer, parameter :: RSEED0 = -777
  integer :: mpierr, rank, npes
  integer :: seed_template(NRANDSEED)

  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, npes, mpierr)
  if (rank == 0) write(*,'(a,i0)') '=== compare_pert, MPI ranks: ', npes

  ! Initialize one RNG state and give an identical copy to each case.  The
  ! generators consume different numbers of variates internally, so they
  ! must not share the same mutable state while running.
  seed_template = 0
  seed_template(1) = RSEED0
  call init_randseed(seed_template)

  call run_latlon(MPI_COMM_WORLD, rank, seed_template)
  call run_sphere(MPI_COMM_WORLD, rank, seed_template)
  call MPI_Finalize(mpierr)

contains

  subroutine run_latlon(comm, rank, seed_template)
    integer, intent(in) :: comm, rank
    integer, intent(in) :: seed_template(NRANDSEED)
    integer :: mpierr, unit
    integer :: rseed(NRANDSEED)
    real :: rfield(N_LON,N_LAT), rfield2(N_LON,N_LAT)
    real(c_double) :: t0, t_init, t_draw
    type(rectangle_random_fields) :: rf
    type(pert_param_type) :: pp
    type(grid_def_type) :: grid

    pp%xcorr = XCORR
    pp%ycorr = YCORR
    pp%coarsen = .false.
    grid%N_lon = N_LON; grid%N_lat = N_LAT
    grid%dlon = DLON; grid%dlat = DLAT
    rseed = seed_template

    call MPI_Barrier(comm, mpierr); t0 = MPI_Wtime()
#ifdef MKL_AVAILABLE
    rf = rectangle_random_fields(pp, grid, comm=comm)
#else
    rf = rectangle_random_fields(pp, grid)
#endif
    call MPI_Barrier(comm, mpierr); t_init = MPI_Wtime() - t0
    call MPI_Barrier(comm, mpierr); t0 = MPI_Wtime()
    call rf%generate_2d_Random_field(rseed, rfield, rfield2, XCORR, YCORR, DLON, DLAT)
    call MPI_Barrier(comm, mpierr); t_draw = MPI_Wtime() - t0

    if (rank == 0) then
      open(newunit=unit, file='latlon_pert.bin', status='replace', access='stream', form='unformatted')
      write(unit) rfield; close(unit)
      call print_stats('latlon', rfield)
      write(*,'(a,f9.4,a,f9.4)') '  init=', t_init, ' draw=', t_draw
    end if
    call rf%finalize()
  end subroutine run_latlon

  subroutine run_sphere(comm, rank, seed_template)
    integer, intent(in) :: comm, rank
    integer, intent(in) :: seed_template(NRANDSEED)
    integer :: mpierr, unit
    integer :: rseed(NRANDSEED)
    real :: rfield(N_LON,N_LAT), rfield2(N_LON,N_LAT)
    real(c_double) :: t0, t_build, t_draw
    type(sphere_random_fields) :: sf
    type(pert_param_type) :: pp
    type(grid_def_type) :: grid

    pp%xcorr = XCORR
    pp%ycorr = YCORR
    grid%N_lon = N_LON; grid%N_lat = N_LAT
    grid%dlon = DLON; grid%dlat = DLAT
    rseed = seed_template

    call MPI_Barrier(comm, mpierr); t0 = MPI_Wtime()
    sf = sphere_random_fields(pp, grid, comm=comm)
    call MPI_Barrier(comm, mpierr); t_build = MPI_Wtime() - t0
    call MPI_Barrier(comm, mpierr); t0 = MPI_Wtime()
    call sf%generate_2d_Random_field(rseed, rfield, rfield2, XCORR, YCORR, DLON, DLAT)
    call MPI_Barrier(comm, mpierr); t_draw = MPI_Wtime() - t0

    if (rank == 0) then
      open(newunit=unit, file='sphere_pert.bin', status='replace', access='stream', form='unformatted')
      write(unit) rfield; close(unit)
      call print_stats('sphere', rfield)
      write(*,'(a,f9.4,a,f9.4)') '  build=', t_build, ' draw=', t_draw
    end if
    call sf%finalize()
  end subroutine run_sphere

  subroutine print_stats(tag, field)
    character(*), intent(in) :: tag
    real, intent(in) :: field(N_LON,N_LAT)
    real :: mean, std
    mean = sum(field) / real(N_LON*N_LAT)
    std = sqrt(sum((field-mean)**2) / real(N_LON*N_LAT-1))
    write(*,'(a,a,a,f9.5,a,f9.5,a,f9.5,a,f9.5)') '  [',tag,'] mean=',mean, &
      ' std=',std,' min=',minval(field),' max=',maxval(field)
  end subroutine print_stats

end program compare_pert
