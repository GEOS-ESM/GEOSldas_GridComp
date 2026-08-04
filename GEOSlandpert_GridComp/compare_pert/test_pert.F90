program test_pert
  use mpi
  use, intrinsic :: iso_fortran_env, only: real32
  use nr_ran2_gasdev, only: NRANDSEED
  use abstract_random_fieldsMod, only: abstract_random_fields
  use rectangle_random_fieldsMod, only: rectangle_random_fields
  use sphere_random_fields_mod, only: sphere_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  implicit none

  integer, parameter :: N_LON = 1440
  integer, parameter :: N_LAT = 720
  integer, parameter :: N_SAMPLES_RECT   = 32
  integer, parameter :: N_SAMPLES_SPHERE = 64
  integer, parameter :: EQUATOR_BAND = 2
  integer, parameter :: RSEED0 = -1701
  real, parameter :: XCORR = 1.0
  real, parameter :: XCORR_SPHERE = 10.0
  real, parameter :: XCORR_TOL        = 0.5   ! tolerance for rectangle (flat FFT)
  real, parameter :: XCORR_TOL_SPHERE = 1.5   ! SMERFS calibration bias ~0.9 deg at 10-deg scale
  real, parameter :: DLON = 360.0 / real(N_LON)
  real, parameter :: DLAT = 180.0 / real(N_LAT)
  real, parameter :: E_FOLD = exp(-0.5)

  integer :: comm, mpierr, rank, npes
  logical :: passed_rectangle, passed_sphere
  type(pert_param_type) :: pp, pp_sphere
  type(grid_def_type) :: grid
  type(rectangle_random_fields) :: rf
  type(sphere_random_fields) :: sf

  comm = MPI_COMM_WORLD
  call MPI_Init(mpierr)
  call MPI_Comm_rank(comm, rank, mpierr)
  call MPI_Comm_size(comm, npes, mpierr)

  pp%xcorr = XCORR
  pp%ycorr = XCORR
  pp%coarsen = .false.
  grid%N_lon = N_LON
  grid%N_lat = N_LAT
  grid%dlon = DLON
  grid%dlat = DLAT

#ifdef MKL_AVAILABLE
  rf = rectangle_random_fields(pp, grid, comm=comm)
#else
  rf = rectangle_random_fields(pp, grid)
#endif
  pp_sphere = pp
  pp_sphere%xcorr = XCORR_SPHERE
  pp_sphere%ycorr = XCORR_SPHERE
  sf = sphere_random_fields(pp_sphere, grid, comm=comm)

  call run_case(rf, 'rectangle', 0, N_SAMPLES_RECT, XCORR, XCORR_TOL, N_LON, N_LAT, &
       DLON, DLAT, passed_rectangle)
  call run_case(sf, 'sphere SMERFS', 50000, N_SAMPLES_SPHERE, XCORR_SPHERE, XCORR_TOL_SPHERE, &
       N_LON, N_LAT, DLON, DLAT, passed_sphere)

  call rf%finalize()
  call sf%finalize()
  call MPI_Finalize(mpierr)

  if (rank == 0 .and. (.not. passed_rectangle .or. .not. passed_sphere)) error stop 1

contains

  subroutine run_case(generator, label, seed_offset, nsamples, target_xcorr, &
       xcorr_tol, nlon, nlat, dlon, dlat, passed_case)
    class(abstract_random_fields), intent(inout) :: generator
    character(*), intent(in) :: label
    integer, intent(in) :: seed_offset, nsamples, nlon, nlat
    real, intent(in) :: target_xcorr, xcorr_tol
    real, intent(in) :: dlon, dlat
    logical, intent(out) :: passed_case
    integer :: sample, lag, i, j, max_lag
    integer :: j_first, j_last, seed_value
    integer :: rseed(NRANDSEED)
    integer :: local_mpierr
    real, allocatable :: rfield(:,:), rfield2(:,:)
    real(real32), allocatable :: corr_sum(:), corr_global(:)
    real :: corr, corr_prev, distance, distance_prev, estimated_xcorr
    logical :: found_crossing

    max_lag = min(nlon / 4, max(1, nint(3.0 * target_xcorr / dlon)))
    allocate(corr_sum(0:max_lag), corr_global(0:max_lag))
    allocate(rfield(nlon, nlat), rfield2(nlon, nlat))
    corr_sum = 0.0_real32
    passed_case = .false.
    j_first = max(1, nlat / 2 - EQUATOR_BAND)
    j_last = min(nlat, nlat / 2 + EQUATOR_BAND)

    do sample = 1, nsamples
       rseed = 0
       seed_value = RSEED0 + seed_offset - rank - (sample - 1) * npes
       rseed(1) = seed_value
       call generator%generate_2d_Random_field(rseed, rfield, rfield2, &
            target_xcorr, target_xcorr, dlon, dlat)
       do lag = 0, max_lag
          do j = j_first, j_last
            do i = 1, nlon
               corr_sum(lag) = corr_sum(lag) + real( &
                     rfield(i,j) * rfield(modulo(i - 1 + lag, nlon) + 1,j), real32)
             end do
          end do
       end do
    end do

    call MPI_Reduce(corr_sum, corr_global, max_lag + 1, MPI_REAL, MPI_SUM, &
         0, comm, local_mpierr)

    if (rank == 0) then
       corr_global = corr_global / corr_global(0)
       found_crossing = .false.
       estimated_xcorr = -1.0
       corr_prev = corr_global(0)
       distance_prev = 0.0

       write(*,'(a,a,a)') '=== test_pert: ', trim(label), ' spatial correlation ==='
       write(*,'(a,f8.3,a,i0,a,i0)') '  target xcorr=', target_xcorr, ' deg, samples=', &
            nsamples * npes, ', equator rows=', j_last - j_first + 1
       write(*,'(a)') '  lag(deg)   correlation'
       do lag = 0, max_lag
          distance = real(lag) * dlon
          corr = corr_global(lag)
          write(*,'(f9.3, f14.6)') distance, corr
          if (.not. found_crossing .and. lag > 0 .and. corr <= E_FOLD) then
             if (corr_prev /= E_FOLD) then
                estimated_xcorr = distance_prev + (E_FOLD - corr_prev) * &
                     (distance - distance_prev) / (corr - corr_prev)
             else
                estimated_xcorr = distance_prev
             end if
             found_crossing = .true.
          end if
          corr_prev = corr
          distance_prev = distance
       end do

       passed_case = found_crossing .and. &
            abs(estimated_xcorr - target_xcorr) <= xcorr_tol
       if (found_crossing) then
          write(*,'(a,f8.3,a,f8.3,a)') '  estimated xcorr=', estimated_xcorr, &
               ' deg (error=', estimated_xcorr - target_xcorr, ' deg)'
       else
          write(*,'(a)') '  estimated xcorr: not found in tested lag range'
       end if
       if (passed_case) then
          write(*,'(a)') '  RESULT: PASS'
       else
          write(*,'(a)') '  RESULT: FAIL'
       end if
    end if
    deallocate(rfield, rfield2, corr_sum, corr_global)
  end subroutine run_case
end program test_pert
