program pert_statistics
  use mpi
  use, intrinsic :: iso_fortran_env, only: real32
  use nr_ran2_gasdev, only: NRANDSEED
  use abstract_random_fieldsMod, only: abstract_random_fields
  use rectangle_random_fieldsMod, only: rectangle_random_fields
  use sphere_random_fields_mod, only: sphere_random_fields
  use LDAS_PertTypes, only: pert_param_type
  use LDAS_TileCoordType, only: grid_def_type
  implicit none

  integer, parameter :: N_LON = 720
  integer, parameter :: N_LAT = 360
  integer, parameter :: N_SAMPLES = 64
  integer, parameter :: EQUATOR_BAND = 2
  integer, parameter :: RSEED0 = -1701
  real, parameter :: TARGET_XCORR(2) = [2.0, 10.0]
  real, parameter :: DLON = 360.0 / real(N_LON)
  real, parameter :: DLAT = 180.0 / real(N_LAT)
  real, parameter :: E_FOLD = exp(-0.5)

  integer :: comm, mpierr, rank, npes, icase
  type(pert_param_type) :: pp
  type(grid_def_type) :: grid
  type(rectangle_random_fields) :: rf
  type(sphere_random_fields) :: sf

  comm = MPI_COMM_WORLD
  call MPI_Init(mpierr)
  call MPI_Comm_rank(comm, rank, mpierr)
  call MPI_Comm_size(comm, npes, mpierr)

  grid%N_lon = N_LON
  grid%N_lat = N_LAT
  grid%dlon = DLON
  grid%dlat = DLAT

  if (rank == 0) then
     write(*,'(a)') '=== pert_statistics: recovered xcorr from generated fields ==='
     write(*,'(a,i0,a,i0,a,i0)') '  grid=', N_LON, 'x', N_LAT, ', MPI ranks=', npes
     write(*,'(a,i0)') '  samples per rank=', N_SAMPLES
  end if

  do icase = 1, size(TARGET_XCORR)
     pp%xcorr = TARGET_XCORR(icase)
     pp%ycorr = TARGET_XCORR(icase)
     pp%coarsen = .false.

#ifdef MKL_AVAILABLE
     rf = rectangle_random_fields(pp, grid, comm=comm)
#else
     rf = rectangle_random_fields(pp, grid)
#endif
     call estimate_xcorr(rf, 'rectangle_random_fields', icase * 10000, TARGET_XCORR(icase))
     call rf%finalize()

     sf = sphere_random_fields(pp, grid, comm=comm)
     call estimate_xcorr(sf, 'sphere_random_fields', 500000 + icase * 10000, TARGET_XCORR(icase))
     call sf%finalize()
  end do

  call MPI_Finalize(mpierr)

contains

  subroutine estimate_xcorr(generator, label, seed_offset, target_xcorr)
    class(abstract_random_fields), intent(inout) :: generator
    character(*), intent(in) :: label
    integer, intent(in) :: seed_offset
    real, intent(in) :: target_xcorr
    integer :: sample, lag, i, j, max_lag
    integer :: j_first, j_last, seed_value
    integer :: rseed(NRANDSEED)
    integer :: local_mpierr
    real, allocatable :: rfield(:,:), rfield2(:,:)
    real(real32), allocatable :: corr_sum(:), corr_global(:)
    real :: corr, corr_prev, distance, distance_prev, estimated_xcorr
    logical :: found_crossing

    max_lag = min(N_LON / 4, max(1, nint(4.0 * target_xcorr / DLON)))
    allocate(corr_sum(0:max_lag), corr_global(0:max_lag))
    allocate(rfield(N_LON, N_LAT), rfield2(N_LON, N_LAT))
    corr_sum = 0.0_real32
    j_first = max(1, N_LAT / 2 - EQUATOR_BAND)
    j_last = min(N_LAT, N_LAT / 2 + EQUATOR_BAND)

    do sample = 1, N_SAMPLES
       rseed = 0
       seed_value = RSEED0 + seed_offset - rank - (sample - 1) * npes
       rseed(1) = seed_value
       call generator%generate_2d_Random_field(rseed, rfield, rfield2, &
            target_xcorr, target_xcorr, DLON, DLAT)
       do lag = 0, max_lag
          do j = j_first, j_last
             do i = 1, N_LON
                corr_sum(lag) = corr_sum(lag) + real( &
                     rfield(i,j) * rfield(modulo(i - 1 + lag, N_LON) + 1,j), real32)
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

       do lag = 0, max_lag
          distance = real(lag) * DLON
          corr = corr_global(lag)
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

       write(*,'(a,a,a,f6.2,a,i0,a)') '  ', trim(label), ': target=', &
            target_xcorr, ' deg, samples=', N_SAMPLES * npes, ','
       if (found_crossing) then
          write(*,'(a,f8.3,a,f8.3,a)') '    recovered xcorr=', estimated_xcorr, &
               ' deg (error=', estimated_xcorr - target_xcorr, ' deg)'
       else
          write(*,'(a)') '    recovered xcorr: not found in tested lag range'
       end if
    end if

    deallocate(rfield, rfield2, corr_sum, corr_global)
  end subroutine estimate_xcorr

end program pert_statistics
