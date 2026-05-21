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

  use mwRTM_routines,                   ONLY:     &
       mwRTM_get_lr_reflectivity

  implicit none

  private

  public :: cygnss_preproc_get_obs_pred

  logical                  :: is_loaded = .false.
  character(len=400)       :: loaded_fname = ''

  integer                  :: N_obs_file = 0
  integer                  :: N_support_file = 0

  integer, allocatable     :: obs_tilenum(:)
  integer, allocatable     :: tile_start(:)
  integer, allocatable     :: tile_count(:)
  integer, allocatable     :: support_tile_index0(:)

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

  subroutine cygnss_preproc_make_fname(this_obs_param, fname)

    ! Build the coefficient-product filename from obs_param%path and
    ! obs_param%name.

    implicit none

    type(obs_param_type), intent(in)  :: this_obs_param
    character(len=400),  intent(out) :: fname

    integer :: lpath

    character(len=*), parameter :: Iam = 'cygnss_preproc_make_fname'

    fname = ''

    if ((len_trim(this_obs_param%path) == 0) .or. (len_trim(this_obs_param%name) == 0)) then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'empty CYGNSS coefficient path/name')
    end if

    lpath = len_trim(this_obs_param%path)

    if (this_obs_param%path(lpath:lpath) == '/') then
       fname = trim(this_obs_param%path) // trim(this_obs_param%name)
    else
       fname = trim(this_obs_param%path) // '/' // trim(this_obs_param%name)
    end if

  end subroutine cygnss_preproc_make_fname

  ! *****************************************************************

  subroutine cygnss_preproc_clear()

    ! Clear the one-file coefficient-product cache.

    implicit none

    if (allocated(obs_tilenum))                   deallocate(obs_tilenum)
    if (allocated(tile_start))                    deallocate(tile_start)
    if (allocated(tile_count))                    deallocate(tile_count)
    if (allocated(support_tile_index0))           deallocate(support_tile_index0)
    if (allocated(coefficient))                   deallocate(coefficient)
    if (allocated(sp_inc_angle))                  deallocate(sp_inc_angle)
    if (allocated(sp_nearest_tile_distance_km))   deallocate(sp_nearest_tile_distance_km)

    is_loaded = .false.
    loaded_fname = ''
    N_obs_file = 0
    N_support_file = 0

  end subroutine cygnss_preproc_clear

  ! *****************************************************************

  subroutine cygnss_preproc_load(this_obs_param)

    ! Load the coefficient-product metadata and sparse support arrays.
    !
    ! The cache intentionally holds one product file at a time.  Repeated
    ! calls with the same obs_param%path/name are cheap.

    implicit none

    type(obs_param_type), intent(in) :: this_obs_param

    integer :: ncid, dimid, varid
    integer :: status

    integer, allocatable :: tmp_tilenum0(:)

    character(len=400) :: fname
    character(len=256) :: product_type
    character(len=256) :: schema_version

    character(len=*), parameter :: Iam = 'cygnss_preproc_load'

    call cygnss_preproc_make_fname(this_obs_param, fname)

    if (is_loaded .and. (trim(fname) == trim(loaded_fname))) return

    call cygnss_preproc_clear()

    if (logit) write(logunit,'(400A)') 'Reading CYGNSS coefficient operator file: ' // trim(fname)

    status = nf90_open(trim(fname), nf90_nowrite, ncid)
    call cygnss_preproc_nc_check(status, Iam, 'opening ' // trim(fname))

    status = nf90_get_att(ncid, nf90_global, 'product_type', product_type)
    call cygnss_preproc_nc_check(status, Iam, 'reading global attribute product_type')

    status = nf90_get_att(ncid, nf90_global, 'schema_version', schema_version)
    call cygnss_preproc_nc_check(status, Iam, 'reading global attribute schema_version')

    if (trim(product_type) /= 'cygnss_tile_coefficient_preprocessor_netcdf') then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unexpected CYGNSS coefficient product_type')
    end if

    if (trim(schema_version) /= '0.3') then
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unexpected CYGNSS coefficient schema_version')
    end if

    status = nf90_inq_dimid(ncid, 'obs', dimid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring obs dimension')

    status = nf90_inquire_dimension(ncid, dimid, len=N_obs_file)
    call cygnss_preproc_nc_check(status, Iam, 'reading obs dimension')

    status = nf90_inq_dimid(ncid, 'support', dimid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring support dimension')

    status = nf90_inquire_dimension(ncid, dimid, len=N_support_file)
    call cygnss_preproc_nc_check(status, Iam, 'reading support dimension')

    allocate(tmp_tilenum0(N_obs_file))
    allocate(obs_tilenum(N_obs_file))
    allocate(tile_start(N_obs_file))
    allocate(tile_count(N_obs_file))
    allocate(sp_inc_angle(N_obs_file))
    allocate(sp_nearest_tile_distance_km(N_obs_file))
    allocate(support_tile_index0(N_support_file))
    allocate(coefficient(N_support_file))

    status = nf90_inq_varid(ncid, 'sp_nearest_tile_index0', varid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring sp_nearest_tile_index0')
    status = nf90_get_var(ncid, varid, tmp_tilenum0)
    call cygnss_preproc_nc_check(status, Iam, 'reading sp_nearest_tile_index0')

    ! The product stores zero-based full-domain tile indices; GEOSldas
    ! observation tile numbers are one-based.

    obs_tilenum = tmp_tilenum0 + 1

    status = nf90_inq_varid(ncid, 'tile_start', varid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring tile_start')
    status = nf90_get_var(ncid, varid, tile_start)
    call cygnss_preproc_nc_check(status, Iam, 'reading tile_start')

    status = nf90_inq_varid(ncid, 'tile_count', varid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring tile_count')
    status = nf90_get_var(ncid, varid, tile_count)
    call cygnss_preproc_nc_check(status, Iam, 'reading tile_count')

    status = nf90_inq_varid(ncid, 'sp_inc_angle', varid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring sp_inc_angle')
    status = nf90_get_var(ncid, varid, sp_inc_angle)
    call cygnss_preproc_nc_check(status, Iam, 'reading sp_inc_angle')

    status = nf90_inq_varid(ncid, 'sp_nearest_tile_distance_km', varid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring sp_nearest_tile_distance_km')
    status = nf90_get_var(ncid, varid, sp_nearest_tile_distance_km)
    call cygnss_preproc_nc_check(status, Iam, 'reading sp_nearest_tile_distance_km')

    status = nf90_inq_varid(ncid, 'tile_index0', varid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring tile_index0')
    status = nf90_get_var(ncid, varid, support_tile_index0)
    call cygnss_preproc_nc_check(status, Iam, 'reading tile_index0')

    status = nf90_inq_varid(ncid, 'coefficient', varid)
    call cygnss_preproc_nc_check(status, Iam, 'inquiring coefficient')
    status = nf90_get_var(ncid, varid, coefficient)
    call cygnss_preproc_nc_check(status, Iam, 'reading coefficient')

    status = nf90_close(ncid)
    call cygnss_preproc_nc_check(status, Iam, 'closing ' // trim(fname))

    deallocate(tmp_tilenum0)

    loaded_fname = fname
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

  integer function cygnss_preproc_find_tile(N_catlH, tile_coord_lH, tilenum)

    ! Map a one-based full-domain tile number to its local-plus-halo index.
    ! The local-plus-halo tile coordinates keep that number in %f_num.

    implicit none

    integer,                                  intent(in) :: N_catlH
    type(tile_coord_type), dimension(N_catlH), intent(in) :: tile_coord_lH
    integer,                                  intent(in) :: tilenum

    integer :: i

    cygnss_preproc_find_tile = -1

    do i=1,N_catlH

       if (tile_coord_lH(i)%f_num == tilenum) then
          cygnss_preproc_find_tile = i
          return
       end if

    end do

  end function cygnss_preproc_find_tile

  ! *****************************************************************

  subroutine cygnss_preproc_get_obs_pred(                                      &
       this_obs_param, N_catlH, tile_coord_lH, N_ens,                          &
       sfmc_lH, mwp_clay_lH, mwp_poros_lH, tilenum, obs_pred)

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
    real, dimension(N_ens),                  intent(out) :: obs_pred

    integer :: obs_ind
    integer :: n_e, k, k1, k2
    integer :: support_tilenum, ind_lH

    real    :: refl_lr, hx_linear

    character(len=400) :: err_msg

    character(len=*), parameter :: Iam = 'cygnss_preproc_get_obs_pred'

    call cygnss_preproc_load(this_obs_param)

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

          ! Support tile indices are zero-based in the coefficient product.
          ! Abort during development if the configured halo does not contain
          ! an explicit support tile needed by this observation.

          support_tilenum = support_tile_index0(k) + 1
          ind_lH = cygnss_preproc_find_tile(N_catlH, tile_coord_lH, support_tilenum)

          if (ind_lH < 1) then
             write(err_msg,*) 'CYGNSS support tile not in halo, owner=', tilenum, &
                  ' support=', support_tilenum
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
