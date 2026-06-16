module cygnss_preprocessed_obs

  ! Read and cache a preprocessed CYGNSS coefficient product, then evaluate
  ! the scalar observation operator:
  !
  !     H(x) = sum_t C_t * R_t(x)
  !
  ! where C_t comes from the coefficient product and R_t is the model LR
  ! reflectivity for each support tile.
  !
  ! amfox+codex, 20 May 2026

  use netcdf

  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_PI

  use enkf_types,                       ONLY:     &
       obs_param_type

  use LDAS_TilecoordType,               ONLY:     &
       tile_coord_type

  use LDAS_ensdrv_globals,              ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_generic,                            &
       LDAS_is_nodata

  use LDAS_exceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use LDAS_DateTimeMod,                 ONLY:     &
       date_time_type,                            &
       augment_date_time,                         &
       datetime_le_refdatetime

  use mwRTM_routines,                   ONLY:     &
       mwRTM_get_lr_reflectivity

  implicit none

  private

  public :: cygnss_preproc_get_obs_pred

  logical                  :: is_loaded = .false.
  character(len=2000)      :: loaded_fname = ''

  integer                  :: N_obs_file = 0
  integer                  :: N_support_file = 0

  integer, allocatable     :: obs_tilenum(:)
  integer, allocatable     :: tile_start(:)
  integer, allocatable     :: tile_count(:)
  integer, allocatable     :: support_tile_ig(:)
  integer, allocatable     :: support_tile_jg(:)

  real,    allocatable     :: coefficient(:)
  real,    allocatable     :: sp_inc_angle(:)
  real,    allocatable     :: sp_nearest_tile_distance_km(:)

contains

  ! *****************************************************************

  subroutine cygnss_preproc_nc_check(status, caller, context)

    ! Check NetCDF return codes and stop with enough context to identify
    ! which product read failed.

    implicit none

    integer,      intent(in) :: status
    character(*), intent(in) :: caller
    character(*), intent(in) :: context

    character(len=400) :: err_msg

    if (status /= nf90_noerr) then
       err_msg = trim(context) // ': ' // trim(nf90_strerror(status))
       call ldas_abort(LDAS_GENERIC_ERROR, caller, err_msg)
    end if

  end subroutine cygnss_preproc_nc_check

  ! *****************************************************************

  subroutine cygnss_preproc_replace_token(string, token, value)

    ! Replace all occurrences of a fixed-width date token in string.
    ! Keep synchronized with read_obs_cygnss_l1_replace_token().

    implicit none

    character(len=*), intent(inout) :: string
    character(len=*), intent(in)    :: token
    character(len=*), intent(in)    :: value

    integer :: ind

    ind = index(string, token)

    do while (ind > 0)
       string = string(:ind-1) // value // string(ind+len(token):)
       ind = index(string, token)
    end do

  end subroutine cygnss_preproc_replace_token

  ! *****************************************************************

  subroutine cygnss_preproc_make_fname(this_obs_param, date_time_file, fname)

    ! Build the coefficient-product filename from obs_param%path and
    ! obs_param%name.  Daily products live under Yyyyy/Mmm directories.

    implicit none

    type(obs_param_type), intent(in)      :: this_obs_param
    type(date_time_type), intent(in)      :: date_time_file
    character(len=400),   intent(out)     :: fname

    integer :: lpath

    character(len=4) :: yyyy
    character(len=2) :: mm, dd
    character(len=400) :: filename

    character(len=*), parameter :: Iam = 'cygnss_preproc_make_fname'

    fname = ''

    if ((len_trim(this_obs_param%path) == 0) .or. (len_trim(this_obs_param%name) == 0)) then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'empty CYGNSS coefficient path/name')
    end if

    lpath = len_trim(this_obs_param%path)
    write(yyyy,'(I4.4)') date_time_file%year
    write(mm,  '(I2.2)') date_time_file%month
    write(dd,  '(I2.2)') date_time_file%day

    filename = trim(this_obs_param%name)
    call cygnss_preproc_replace_token(filename, 'yyyymmdd', yyyy//mm//dd)

    if (this_obs_param%path(lpath:lpath) == '/') then
       fname = trim(this_obs_param%path) // 'Y' // yyyy // '/M' // mm // '/' // trim(filename)
    else
       fname = trim(this_obs_param%path) // '/Y' // yyyy // '/M' // mm // '/' // trim(filename)
    end if

  end subroutine cygnss_preproc_make_fname

  ! *****************************************************************

  subroutine cygnss_preproc_clear()

    ! Clear the daily-file-set coefficient-product cache.

    implicit none

    if (allocated(obs_tilenum))                   deallocate(obs_tilenum)
    if (allocated(tile_start))                    deallocate(tile_start)
    if (allocated(tile_count))                    deallocate(tile_count)
    if (allocated(support_tile_ig))               deallocate(support_tile_ig)
    if (allocated(support_tile_jg))               deallocate(support_tile_jg)
    if (allocated(coefficient))                   deallocate(coefficient)
    if (allocated(sp_inc_angle))                  deallocate(sp_inc_angle)
    if (allocated(sp_nearest_tile_distance_km))   deallocate(sp_nearest_tile_distance_km)

    is_loaded = .false.
    loaded_fname = ''
    N_obs_file = 0
    N_support_file = 0

  end subroutine cygnss_preproc_clear

  ! *****************************************************************

  subroutine cygnss_preproc_load(this_obs_param, date_time, dtstep_assim)

    ! Load the coefficient-product metadata and sparse support arrays.
    !
    ! The cache holds all daily product files touched by the current
    ! assimilation window.  Repeated calls for the same window are cheap.

    implicit none

    type(obs_param_type), intent(in)    :: this_obs_param
    type(date_time_type), intent(in)    :: date_time
    integer,             intent(in)    :: dtstep_assim

    integer :: ncid, dimid, varid
    integer :: status
    integer :: N_obs_this, N_support_this
    integer :: N_obs_total, N_support_total
    integer :: obs_offset, support_offset

    integer, allocatable :: tmp_tilenum0(:)
    integer, allocatable :: tmp_tile_start(:)

    character(len=400) :: fname
    character(len=2000) :: fname_key
    character(len=256) :: product_type
    character(len=256) :: schema_version

    logical :: file_exists
    type(date_time_type) :: date_time_low, date_time_up, date_time_file

    character(len=*), parameter :: Iam = 'cygnss_preproc_load'

    date_time_low = date_time
    call augment_date_time( -(dtstep_assim/2), date_time_low )
    date_time_up = date_time
    call augment_date_time(  (dtstep_assim/2), date_time_up )

    date_time_file = date_time_type(date_time_low%year, date_time_low%month, date_time_low%day, &
         0, 0, 0, -9999, -9999)

    fname_key = ''

    ! Cheap cache-key pass.  Avoid opening NetCDF files on every warm-cache
    ! observation lookup; the expensive dimension-count pass is only needed
    ! when the date-window file set changes.
    do while (datetime_le_refdatetime(date_time_file, date_time_up))

       call cygnss_preproc_make_fname(this_obs_param, date_time_file, fname)
       inquire(file=trim(fname), exist=file_exists)

       if (file_exists) then
          fname_key = trim(fname_key) // '|' // trim(fname)
       end if

       call augment_date_time(86400, date_time_file)

    end do

    if (is_loaded .and. (trim(fname_key) == trim(loaded_fname))) return

    call cygnss_preproc_clear()

    N_obs_total = 0
    N_support_total = 0

    date_time_file = date_time_type(date_time_low%year, date_time_low%month, date_time_low%day, &
         0, 0, 0, -9999, -9999)

    do while (datetime_le_refdatetime(date_time_file, date_time_up))

       call cygnss_preproc_make_fname(this_obs_param, date_time_file, fname)
       inquire(file=trim(fname), exist=file_exists)

       if (file_exists) then

          status = nf90_open(trim(fname), nf90_nowrite, ncid)
          call cygnss_preproc_nc_check(status, Iam, 'opening ' // trim(fname))

          status = nf90_get_att(ncid, nf90_global, 'product_type', product_type)
          call cygnss_preproc_nc_check(status, Iam, 'reading global attribute product_type')

          status = nf90_get_att(ncid, nf90_global, 'schema_version', schema_version)
          call cygnss_preproc_nc_check(status, Iam, 'reading global attribute schema_version')

          if (trim(product_type) /= 'cygnss_tile_coefficient_preprocessor_netcdf') then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unexpected CYGNSS coefficient product_type')
          end if

          if (trim(schema_version) /= '0.4') then
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, &
                  'CYGNSS coefficient schema_version must be 0.4 - regenerate files with preprocessor >= 0.4')
          end if

          status = nf90_inq_dimid(ncid, 'obs', dimid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring obs dimension')
          status = nf90_inquire_dimension(ncid, dimid, len=N_obs_this)
          call cygnss_preproc_nc_check(status, Iam, 'reading obs dimension')

          status = nf90_inq_dimid(ncid, 'support', dimid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring support dimension')
          status = nf90_inquire_dimension(ncid, dimid, len=N_support_this)
          call cygnss_preproc_nc_check(status, Iam, 'reading support dimension')

          status = nf90_close(ncid)
          call cygnss_preproc_nc_check(status, Iam, 'closing ' // trim(fname))

          N_obs_total = N_obs_total + N_obs_this
          N_support_total = N_support_total + N_support_this

       end if

       call augment_date_time(86400, date_time_file)

    end do

    if (N_obs_total < 1) then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'no CYGNSS coefficient operator files found for window')
    end if

    N_obs_file = N_obs_total
    N_support_file = N_support_total

    allocate(obs_tilenum(N_obs_file))
    allocate(tile_start(N_obs_file))
    allocate(tile_count(N_obs_file))
    allocate(sp_inc_angle(N_obs_file))
    allocate(sp_nearest_tile_distance_km(N_obs_file))
    allocate(support_tile_ig(N_support_file))
    allocate(support_tile_jg(N_support_file))
    allocate(coefficient(N_support_file))

    obs_offset = 0
    support_offset = 0

    date_time_file = date_time_type(date_time_low%year, date_time_low%month, date_time_low%day, &
         0, 0, 0, -9999, -9999)

    do while (datetime_le_refdatetime(date_time_file, date_time_up))

       call cygnss_preproc_make_fname(this_obs_param, date_time_file, fname)
       inquire(file=trim(fname), exist=file_exists)

       if (file_exists) then

          if (logit) write(logunit,'(400A)') 'Reading CYGNSS coefficient operator file: ' // trim(fname)

          status = nf90_open(trim(fname), nf90_nowrite, ncid)
          call cygnss_preproc_nc_check(status, Iam, 'opening ' // trim(fname))

          status = nf90_inq_dimid(ncid, 'obs', dimid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring obs dimension')
          status = nf90_inquire_dimension(ncid, dimid, len=N_obs_this)
          call cygnss_preproc_nc_check(status, Iam, 'reading obs dimension')

          status = nf90_inq_dimid(ncid, 'support', dimid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring support dimension')
          status = nf90_inquire_dimension(ncid, dimid, len=N_support_this)
          call cygnss_preproc_nc_check(status, Iam, 'reading support dimension')

          allocate(tmp_tilenum0(N_obs_this))
          allocate(tmp_tile_start(N_obs_this))

          status = nf90_inq_varid(ncid, 'sp_nearest_tile_index0', varid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring sp_nearest_tile_index0')
          status = nf90_get_var(ncid, varid, tmp_tilenum0)
          call cygnss_preproc_nc_check(status, Iam, 'reading sp_nearest_tile_index0')
          obs_tilenum(obs_offset+1:obs_offset+N_obs_this) = tmp_tilenum0 + 1

          status = nf90_inq_varid(ncid, 'tile_start', varid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring tile_start')
          status = nf90_get_var(ncid, varid, tmp_tile_start)
          call cygnss_preproc_nc_check(status, Iam, 'reading tile_start')
          tile_start(obs_offset+1:obs_offset+N_obs_this) = tmp_tile_start + support_offset

          status = nf90_inq_varid(ncid, 'tile_count', varid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring tile_count')
          status = nf90_get_var(ncid, varid, tile_count(obs_offset+1:obs_offset+N_obs_this))
          call cygnss_preproc_nc_check(status, Iam, 'reading tile_count')

          status = nf90_inq_varid(ncid, 'sp_inc_angle', varid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring sp_inc_angle')
          status = nf90_get_var(ncid, varid, sp_inc_angle(obs_offset+1:obs_offset+N_obs_this))
          call cygnss_preproc_nc_check(status, Iam, 'reading sp_inc_angle')

          status = nf90_inq_varid(ncid, 'sp_nearest_tile_distance_km', varid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring sp_nearest_tile_distance_km')
          status = nf90_get_var(ncid, varid, sp_nearest_tile_distance_km(obs_offset+1:obs_offset+N_obs_this))
          call cygnss_preproc_nc_check(status, Iam, 'reading sp_nearest_tile_distance_km')

          status = nf90_inq_varid(ncid, 'tile_ig', varid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring tile_ig')
          status = nf90_get_var(ncid, varid, support_tile_ig(support_offset+1:support_offset+N_support_this))
          call cygnss_preproc_nc_check(status, Iam, 'reading tile_ig')

          status = nf90_inq_varid(ncid, 'tile_jg', varid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring tile_jg')
          status = nf90_get_var(ncid, varid, support_tile_jg(support_offset+1:support_offset+N_support_this))
          call cygnss_preproc_nc_check(status, Iam, 'reading tile_jg')

          status = nf90_inq_varid(ncid, 'coefficient', varid)
          call cygnss_preproc_nc_check(status, Iam, 'inquiring coefficient')
          status = nf90_get_var(ncid, varid, coefficient(support_offset+1:support_offset+N_support_this))
          call cygnss_preproc_nc_check(status, Iam, 'reading coefficient')

          status = nf90_close(ncid)
          call cygnss_preproc_nc_check(status, Iam, 'closing ' // trim(fname))

          deallocate(tmp_tilenum0)
          deallocate(tmp_tile_start)

          obs_offset = obs_offset + N_obs_this
          support_offset = support_offset + N_support_this

       end if

       call augment_date_time(86400, date_time_file)

    end do

    loaded_fname = fname_key
    is_loaded = .true.

  end subroutine cygnss_preproc_load

  ! *****************************************************************

  integer function cygnss_preproc_find_obs(tilenum)

    ! Find the cached coefficient-product observation for a GEOSldas owner
    ! tile.  If duplicates exist, use the observation whose specular point is
    ! closest to the owner tile.

    implicit none

    integer, intent(in) :: tilenum

    integer :: i
    real    :: best_distance

    cygnss_preproc_find_obs = -1
    best_distance = huge(best_distance)

    do i=1,N_obs_file

       if (obs_tilenum(i) == tilenum) then

          if (sp_nearest_tile_distance_km(i) < best_distance) then
             best_distance = sp_nearest_tile_distance_km(i)
             cygnss_preproc_find_obs = i
          end if

       end if

    end do

  end function cygnss_preproc_find_obs

  ! *****************************************************************

  integer function cygnss_preproc_find_tile_by_ig_jg(N_catlH, tile_coord_lH, ig, jg)

    ! Map global EASE-grid (ig, jg) coordinates to the local-plus-halo index.
    ! i_indg/j_indg are globally fixed for any M36 tile regardless of experiment.

    implicit none

    integer,                                   intent(in) :: N_catlH
    type(tile_coord_type), dimension(N_catlH), intent(in) :: tile_coord_lH
    integer,                                   intent(in) :: ig, jg

    integer :: i

    cygnss_preproc_find_tile_by_ig_jg = -1

    do i=1,N_catlH

       if (tile_coord_lH(i)%i_indg == ig .and. tile_coord_lH(i)%j_indg == jg) then
          cygnss_preproc_find_tile_by_ig_jg = i
          return
       end if

    end do

  end function cygnss_preproc_find_tile_by_ig_jg

  ! *****************************************************************

  subroutine cygnss_preproc_get_obs_pred(                                      &
       this_obs_param, N_catlH, tile_coord_lH, N_ens,                          &
       sfmc_lH, mwp_clay_lH, mwp_poros_lH, tilenum,                            &
       date_time, dtstep_assim, obs_pred)

    ! Evaluate the CYGNSS coefficient operator for one GEOSldas observation
    ! and all ensemble members.
    !
    ! This operator bypasses the generic FOV weighting in get_obs_pred()
    ! because the preprocessor has already identified the exact support tiles
    ! and coefficient weights for this scalar observation.

    implicit none

    type(obs_param_type),                    intent(in)  :: this_obs_param
    integer,                                 intent(in)  :: N_catlH
    type(tile_coord_type), dimension(N_catlH), intent(in) :: tile_coord_lH
    integer,                                 intent(in)  :: N_ens

    real, dimension(N_catlH,N_ens),          intent(in)  :: sfmc_lH
    real, dimension(N_catlH,N_ens),          intent(in)  :: mwp_clay_lH
    real, dimension(N_catlH,N_ens),          intent(in)  :: mwp_poros_lH

    integer,                                 intent(in)  :: tilenum
    type(date_time_type),                    intent(in)  :: date_time
    integer,                                 intent(in)  :: dtstep_assim
    real, dimension(N_ens),                  intent(out) :: obs_pred

    integer :: obs_ind
    integer :: n_e, k, k1, k2
    integer :: ind_lH

    real    :: refl_lr, hx_linear

    character(len=400) :: err_msg

    character(len=*), parameter :: Iam = 'cygnss_preproc_get_obs_pred'

    call cygnss_preproc_load(this_obs_param, date_time, dtstep_assim)

    obs_pred = nodata_generic

    obs_ind = cygnss_preproc_find_obs(tilenum)

    if (obs_ind < 1) then
       write(err_msg,*) 'CYGNSS coefficient obs not found for tilenum=', tilenum
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    k1 = tile_start(obs_ind) + 1
    k2 = tile_start(obs_ind) + tile_count(obs_ind)

    if ((k1 < 1) .or. (k2 > N_support_file) .or. (k2 < k1)) then
       write(err_msg,*) 'invalid CYGNSS support range for tilenum=', tilenum
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    do n_e=1,N_ens

       hx_linear = 0.

       do k=k1,k2

          ind_lH = cygnss_preproc_find_tile_by_ig_jg(N_catlH, tile_coord_lH, &
               support_tile_ig(k), support_tile_jg(k))

          if (ind_lH < 1) then
             write(err_msg,*) 'CYGNSS support tile not in halo, owner=', tilenum, &
                  ' ig=', support_tile_ig(k), ' jg=', support_tile_jg(k)
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
          end if

          call mwRTM_get_lr_reflectivity(                                       &
               this_obs_param%freq, sp_inc_angle(obs_ind),                      &
               mwp_clay_lH(ind_lH,n_e), mwp_poros_lH(ind_lH,n_e),               &
               sfmc_lH(ind_lH,n_e), refl_lr )

          if (LDAS_is_nodata(refl_lr)) then
             hx_linear = nodata_generic
             exit
          end if

          hx_linear = hx_linear + coefficient(k) * refl_lr

       end do

       if (.not. LDAS_is_nodata(hx_linear)) then
          if (hx_linear > 0.) obs_pred(n_e) = 10. * log10(hx_linear)
       end if

    end do

  end subroutine cygnss_preproc_get_obs_pred

end module cygnss_preprocessed_obs
