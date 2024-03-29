
module mwRTM_routines

  ! subroutines for microwave radiative transfer model
  !
  ! Select a specific configuration of the RTM via the field
  ! "RTM_ID" in the "obs_param" type.
  !
  ! %RTM_ID = ID of radiative transfer model to use for Tb forward modeling
  !           (subroutine get_obs_pred())
  !           0 = none
  !           1 = L-band tau-omega model as in De Lannoy et al. 2013 (doi:10.1175/JHM-D-12-092.1) (SMOS)
  !           2 = same as 1 but without Pellarin atm corr (SMAP)
  !           3 = same as 1 but with Mironov and SMAP L2_SM pol mixing (SMOS)
  !           4 = same as 3 but without Pellarin atm corr (targeted for SMAP L4_SM Version 8)
  !
  ! reichle, 16 May 2011
  ! reichle, 31 Mar 2015 - added RTM_ID
  !
  ! --------------------------------------------------------------------------

  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_PI,                                   &
       MAPL_TICE

  use mwRTM_types,                      ONLY:     &
       mwRTM_param_type,                          &
       mwRTM_param_nodata_check,                  &
       assignment (=)

  use LDAS_ensdrv_globals,              ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_generic,                            &
       nodata_tol_generic,                        &
       LDAS_is_nodata

  use LDAS_exceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  implicit none

  ! everything is private by default unless made public

  private

  public :: mwRTM_get_Tb, catch2mwRTM_vars

  ! ---------------------------------------------------------

  real,    parameter :: Tb_sky      = 2.7        ! cosmic mw background temp  [K]

  complex, parameter :: diel_ice    = (3.2, 0.1) ! dielec. const. of ice
  complex, parameter :: diel_air    = (1.0, 0.0) ! dielec. const. of air
  complex, parameter :: diel_rock   = (5.5, 0.2) ! dielec. const. of rock

  real,    parameter :: diel_watinf =  4.9       ! dielec. const. of water @ high freq.
                                                 ! (Stogryn 1971)

  real,    parameter :: eps_0       = 8.854e-12  ! vacuum permittivity     [Farads/meter]
                                                 ! (Klein and Swift 1977)

  real,    parameter :: rho_soil    = 2.66       ! soil specific density   [g/cm3]

contains

  ! **********************************************************************

  ! Subroutine mwRTM_get_param() reads binary mwRTM files and is no longer used.
  !
  ! The subroutine has been replaced by:
  ! - Components/GEOSldas_GridComp/GEOSldas_App/[..]/mwrtm_bin2nc4.F90 converts mwRTM files from binary to nc4
  ! - get_mwrtm_param() in GEOS_LandAssimGridComp.F90 converts the internal state
  !    variables of the Land Assim GridComp into the mwRTM structure.
  !
  ! reichle, 4 Aug 2020

!  subroutine mwRTM_get_param( N_catg, N_tile, d2g, tile_id, mwRTM_param_path, &
!       need_mwRTM_param, mwp)
!
!    ! Read microwave RTM parameters from file.
!    !
!    ! reichle, 17 May 2011
!    ! reichle, 21 Oct 2011 - added input of mwRTM_param from file
!    ! reichle, 23 Oct 2012 - removed look-up table option (too complicated with
!    !                          "new" (200+) soil classes)
!    !                      - added tile_id check when reading mwRTM params from file
!
!    implicit none
!
!    integer,                                   intent(in)  :: N_catg, N_tile
!
!    integer,                dimension(N_tile), intent(in)  :: d2g, tile_id
!
!    character(200),                            intent(in)  :: mwRTM_param_path
!
!    logical,                                   intent(in)  :: need_mwRTM_param
!
!    type(mwRTM_param_type), dimension(N_tile), intent(out) :: mwp  ! mwRTM parameters
!
!    ! local variables
!
!    integer, parameter                          :: N_search_dir_max = 5
!
!    integer                                     :: n, N_search_dir, istat
!
!    character( 80)                              :: fname
!
!    character(100), dimension(N_search_dir_max) :: search_dir
!
!    logical                                     :: all_nodata, mwp_nodata
!
!    character(len=*), parameter :: Iam = 'mwRTM_get_param'
!    character(len=400) :: err_msg
!
!    ! ----------------------------------------------------------------------
!    !
!    ! initialize
!
!    do n=1,N_tile
!       mwp(n) = nodata_generic
!    end do
!
!    ! read mwRTM parameters from file
!
!    if (logit) write (logunit,*) 'Reading microwave RTM parameters from file'
!
!    fname = '/mwRTM_param.bin'
!
!    N_search_dir = 2  ! specify sub-dirs of mwRTM_param_path to search for file "fname"
!
!    search_dir(1) = 'mwRTM'
!    search_dir(2) = '.'
!
!    ! when called with optional argument "istat" subroutine "open_land_param_file()"
!    ! will *NOT* stop upon failure to open the file
!
!    istat = open_land_param_file( 10, .false., .true., N_search_dir, fname, &
!         mwRTM_param_path, search_dir, ignore_stop=.true.)
!
!    if (istat==0) then
!
!       call io_mwRTM_param_type( 'r', 10, N_tile, mwp, N_catg, tile_id, d2g )
!
!       close (10,status='keep')
!
!       if (logit) write (logunit,*) 'done reading'
!       if (logit) write (logunit,*)
!
!    else
!
!       if (logit) write (logunit,*) 'WARNING: Could not open file!'
!       if (logit) write (logunit,*)
!
!    end if
!
!    ! check for no-data-values in parameters
!    ! if any field is a nodata value, set all fields to nodata value
!
!    all_nodata = .true.
!
!    do n=1,N_tile
!
!       call mwRTM_param_nodata_check( mwp(n), mwp_nodata )
!
!       if (.not. mwp_nodata) all_nodata = .false.
!
!    end do
!
!    ! stop if mwRTM parameters needed but not available (ie, all are no-data)
!
!    if (all_nodata .and. need_mwRTM_param) then
!       err_msg = 'mwRTM params needed but all are no-data!'
!       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
!    end if
!
!    ! warn if mwRTM parameters are all no-data (may not be needed)
!
!    if (all_nodata .and. logit) then
!
!       write (logunit,*) '#########################################################'
!       write (logunit,*)
!       write (logunit,*) '  WARNING: *All* parameters for the microwave radiative  '
!       write (logunit,*) '           transfer model (mwRTM) are no-data values!!!  '
!       write (logunit,*)
!       write (logunit,*) '#########################################################'
!       write (logunit,*)
!
!    end if
!
!  end subroutine mwRTM_get_param

  ! ****************************************************************

  subroutine catch2mwRTM_vars( N_tile, vegcls_catch, poros_catch, poros_mwRTM, &
       sfmc_catch, tsurf_catch, tp1_catch, sfmc_mwRTM, tsoil_mwRTM, tp1_in_Kelvin )

    ! convert soil moisture, surface temperature, and soil temperature from the Catchment
    ! model into soil moisture and soil temperature inputs for the microwave radiative
    ! transfer model (mwRTM)
    !
    ! reichle, 11 Dec 2013
    !
    ! added optional switch to allow tp1_catch input in Kelvin
    ! - reichle & borescan, 6 Nov 2020

    implicit none

    integer,                    intent(in)  :: N_tile

    integer, dimension(N_tile), intent(in)  :: vegcls_catch

    real,    dimension(N_tile), intent(in)  :: poros_catch, poros_mwRTM
    real,    dimension(N_tile), intent(in)  :: sfmc_catch,  tsurf_catch, tp1_catch

    real,    dimension(N_tile), intent(out) :: sfmc_mwRTM,  tsoil_mwRTM

    logical,                    intent(in),  optional :: tp1_in_Kelvin

    ! local variables

    logical                    :: tp1_in_K

    real,    dimension(N_tile) :: tp1

    ! -----------------------------------------------------------------------
    !
    ! reichle, 22 Oct 2012: scaling factor added because it is necessary for
    !                       proper functioning of mwRTM calibration to SMOS obs

    sfmc_mwRTM  = sfmc_catch * poros_mwRTM / poros_catch

    ! diagnose soil temperature to be used with mwRTM
    ! (change prompted by revision of Catchment model parameter CSOIL_2)
    ! - reichle, 23 Dec 2015

    ! NOTE: By default, "tp" is assumed to be in deg Celsius

    tp1_in_K = .false.

    if (present(tp1_in_Kelvin)) tp1_in_K = tp1_in_Kelvin

    tp1 = tp1_catch

    if (.not. tp1_in_K)  tp1 = tp1 + MAPL_TICE  ! convert units if not already in Kelvin

    tsoil_mwRTM = tp1

  end subroutine catch2mwRTM_vars

  ! ****************************************************************

  subroutine mwRTM_get_Tb( N_tile, freq, inc_angle, mwp, elev,             &
       LAI, soilmoist, soiltemp, SWE, Tair, RTM_ID, Tb_h, Tb_v )

    !---------------------------------------------------
    !RTM adapted from Steven Chan
    !23 Nov 2010: - adapted to include the atmospheric correction
    !               GDL, code based on CMEMv3.0
    !
    !Instead of passing on the whole diagnostic structure,
    !only pass rtmv, T5, Tc !GDL, 22Oct10
    !
    !N_cat could be either all tiles or only a part of them!
    !
    !== rtmv: Soil moisture in in g/cm3 (volumetric)
    !FYI:
    !SM vol/vol ~ g/cm3
    !gravim SM (g/g) * bulk densit (g/cm3) = volum SM (vol/vol) ~ g/cm3
    !
    !== VWC: Vegetation water content in kg/m2 (=mm)
    !== T5:  Soil temperature [K]
    !== Tc:  Vegetation canopy temperature [K]
    !
    !21 Mar 2011: - temperature treatment changed
    !==> instead of passing on T5 and Tc (= tsurf)
    !we now pass on tsurf and tp1 and diagnose T10cm (~T5) from it HERE
    !==> tsurf and tp1 are now input variables
    !==> T5 and Tc are local variables
    !
    !    Mar 2011: - replaced realdobson by Wang for diel cst
    ! 02 May 2011: - put in *COS(inc*d2r)**Nrh/v for rsh/v
    !              - changed Q=0 to Q=expression in CMEM
    !              - option to scale model SM before entering RTM
    !
    ! 16 May 2011: - reichle: - included in LDASsa, major revisions
    ! 21 Oct 2011: - reichle: input "poros" now via "mwp"
    ! 23 Nov 2011: - reichle: changed tsoil_threshold b/c QC now done for individual
    !                         ensemble members
    ! 22 Oct 2012: - reichle: removed interception water ("capac") contribution
    !
    !---------------------------------------------------

    implicit none

    integer,                                   intent(in)  :: N_tile     ! number of tiles

    real,                                      intent(in)  :: freq       ! [Hz]
    real,                                      intent(in)  :: inc_angle  ! [deg]

    type(mwRTM_param_type), dimension(N_tile), intent(in)  :: mwp

    real,                   dimension(N_tile), intent(in)  :: elev       ! [m]

    real,                   dimension(N_tile), intent(in)  :: LAI        ! [dim-less]
    real,                   dimension(N_tile), intent(in)  :: soilmoist  ! [m3/m3]
    real,                   dimension(N_tile), intent(in)  :: soiltemp   ! [K]
    real,                   dimension(N_tile), intent(in)  :: SWE        ! [kg/m2] "mm"
    real,                   dimension(N_tile), intent(in)  :: Tair       ! [K]

    integer,                                   intent(in)  :: RTM_ID

    real,                   dimension(N_tile), intent(out) :: Tb_h, Tb_v ! [K]

    ! --------------------

    ! local variables

    real, parameter :: tsoil_threshold = MAPL_TICE+0.2 ! avoid "frozen" soil  [K]

    real, parameter :: SWE_threshold   = 1.e-4         ! avoid snow           [kg/m2]

    integer :: n

    real    :: inc, sin_inc, cos_inc

    complex :: c_er, tmpc1, tmpc2

    real    :: roh, rov, rsh, rsv

    real    :: h_mc, Q, slope

    real    :: vwc, Ah, Av, exptauh, exptauv, exptauh2, exptauv2, tmpreal

    real    :: exptau_atm, tau_atm, Tb_ad, Tb_au

    real    :: soiltemp_in_C, Tc

    !real    ::  er_r   ! for realdobson

    logical :: veg_params_nodata

    character(len=*), parameter :: Iam = 'mwRTM_get_Tb'
    character(len=400) :: err_msg

    !---------------------------------------------------

    !if (logit) write(logunit,*) 'entering mwRTM_get_Tb...'

    ! check first element of elevation against no-data-value
    ! (elevation is only needed for RTMs with Pellarin atm corr)

    select case (RTM_ID)

    case(1,3)

       if ( abs(elev(1)-nodata_generic)<nodata_tol_generic ) then
          err_msg = 'elev is no-data-value'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if

    case(2,4)

       ! do nothing

    case default

       err_msg = 'unknown RTM_ID (during check for elevation)'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end select

    ! pre-compute sine and cosine

    inc = inc_angle * MAPL_PI/180.0   ! degrees to radians

    sin_inc = sin(inc)
    cos_inc = cos(inc)

    do n=1,N_tile

       ! compute Tb only under snow-free and non-frozen conditions
       ! and only where mwRTM parameters are available
       !
       ! abbreviated check for "good" mwRTM parameter values here relies on prep work done
       ! in subroutine mwRTM_param_nodata_check(), which was called in mwRTM_get_param()

       veg_params_nodata = ( LDAS_is_nodata(mwp(n)%bh) .and. LDAS_is_nodata(mwp(n)%vegopacity) )

       if ( (SWE(n)<SWE_threshold)                                .and.         &
            (soiltemp(n)>tsoil_threshold)                         .and.         &
            (.not. LDAS_is_nodata(mwp(n)%sand))                   .and.         &
            (.not. veg_params_nodata)                                           &
            )                                                      then

          ! soil dielectric constant

          soiltemp_in_C = soiltemp(n) - MAPL_TICE

          select case (RTM_ID)

          case(1,2)

             CALL DIELWANG (freq, soilmoist(n),   soiltemp_in_C,        &
                  mwp(n)%wang_wt, mwp(n)%wang_wp, mwp(n)%poros,         &
                  mwp(n)%sand,    mwp(n)%clay,    c_er            )

          case(3,4)

             CALL MIRONOV( freq, soilmoist(n), mwp(n)%clay, c_er )

          case default

             err_msg = 'unknown RTM_ID (during soil dielectric model)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

          end select

          ! soil reflectivity for smooth surface based on dielect const. (Fresnel)

          tmpc1 = SQRT(c_er - sin_inc**2)

          tmpc2 = c_er * cos_inc

          roh   = ABS( (cos_inc - tmpc1) / (cos_inc + tmpc1) )**2
          rov   = ABS( (tmpc2   - tmpc1) / (tmpc2   + tmpc1) )**2

          ! -------------------------------------------------------
          !
          ! roughness corrections:
          !
          ! soil reflectivity for rough surface, based on h-Q model
          ! note that in Choudhury et al., 79, there is a factor exp(-h cos^2 inc)
          ! GDL, 14Feb11, replaced h by SM-dependent h
          ! GDL, 02May11, added cos^N factor
          !
          ! 1) roughness parameter h depends on soil moisture

          if     (soilmoist(n)<=mwp(n)%rgh_wmin) then

             h_mc = mwp(n)%rgh_hmax

          elseif (soilmoist(n)>=mwp(n)%rgh_wmax) then

             h_mc = mwp(n)%rgh_hmin

          else

             slope = &
                  (mwp(n)%rgh_hmin - mwp(n)%rgh_hmax)/ &
                  (mwp(n)%rgh_wmax - mwp(n)%rgh_wmin)

             h_mc = mwp(n)%rgh_hmax + slope * (soilmoist(n) - mwp(n)%rgh_wmin)

          endif

          ! 2) polarization mixing, Q as defined in CMEM:

          if (freq < 2.e9) then

             select case (RTM_ID)

             case(1,2)

                Q = 0.   ! Q is assumed zero at low frequency

             case(3,4)

                Q = 0.1771 * h_mc

             case default

                err_msg = 'unknown RTM_ID (during pol mixing)'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

             end select

          else

             Q = 0.35 * (1.0 - exp(-0.6 * (mwp(n)%rgh_polmix**2) * (freq/1.e9) ))

          end if

          ! rough surface reflectivity

          select case (RTM_ID)

          case(1,2)

             rsh = ( (1-Q) * roh + Q * rov) * EXP(-h_mc*cos_inc**mwp(n)%rgh_nrh)
             rsv = ( (1-Q) * rov + Q * roh) * EXP(-h_mc*cos_inc**mwp(n)%rgh_nrv)

          case(3,4)

             rsh = ( (1-Q) * roh + Q * rov) * EXP(-h_mc*cos_inc**2)
             rsv = ( (1-Q) * rov + Q * roh) * EXP(-h_mc*cos_inc**2)

          case default

             err_msg = 'unknown RTM_ID (during rough reflectivity)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

          end select

          ! -------------------------------------------------------------
          !
          ! Tb at top of vegetation (excl atmos contribution)  (tau-omega model)

          ! vegetation attenuation parameters can come from either of two sources:
          ! - static look-up table or calibrated parameters (bh, bv, lewt), to be combined
          !   with time-varying LAI
          ! - vegetation opacity (from file); varies with time

          ! the if statement above, in conjunction with mwRTM_param_nodata_check() called earlier,
          ! ensures that we have "good" values for either (vegopacity) or (bh,bv,lewt)

          if ( LDAS_is_nodata( mwp(n)%vegopacity ) ) then

             ! == vwc: Vegetation water content in kg/m2 (=mm)
             ! needs to be the total columnar vegetation water content,
             ! depends on greenness/NDVI/LAI...
             ! VWC=LEWT*LAI, lewt is actually a time-varying parameter!
             ! For now LEWT is guessed based on literature, and kept cst.

             vwc = mwp(n)%lewt * lai(n)

             ! removed contribution of interception water ("capac")
             !
             !! ! add bit of intercepted water as well
             !!
             !! vwc = vwc + capac(n)

             ! Vegetation optical thickness tau=b*VWC  (eq. (2) in Crow et al. 2005)

             tmpreal = vwc/cos_inc

             exptauh = EXP( -mwp(n)%bh * tmpreal )
             exptauv = EXP( -mwp(n)%bv * tmpreal )

          else

             exptauh = EXP( -mwp(n)%vegopacity)
             exptauv = EXP( -mwp(n)%vegopacity)

          end if

          Tc = soiltemp(n)        ! canopy temp = soil temp

          tmpreal = Tc * (1. - mwp(n)%omega)

          Ah = tmpreal * (1. - exptauh)
          Av = tmpreal * (1. - exptauv)

          ! Eq.(1) in Crow et al. 2005:
          !
          ! Tb_tov = Tb_soil.e_p.exp(-tau/cos theta) + Tb_veg.(1 + r_r.exp(-tau/cos theta))

          Tb_h(n) = soiltemp(n) * (1. - rsh) * exptauh + Ah * (1. + rsh * exptauh)
          Tb_v(n) = soiltemp(n) * (1. - rsv) * exptauv + Av * (1. + rsv * exptauv)


          ! -------------------------------------------------------------
          !
          ! Atmospheric correction
          !
          ! GDL 23nov10

          select case (RTM_ID)

          case(1,3)

             exptauh2 = exptauh * exptauh
             exptauv2 = exptauv * exptauv

             if (freq<2.e9) then

                call ATMPELLARIN( elev(n)/1000., Tair(n), cos_inc, tau_atm, Tb_ad, Tb_au )

             else

                err_msg = 'cannot compute atm corr for given freq'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

             end if

             exptau_atm = EXP( - tau_atm/cos_inc )

             Tb_h(n) = Tb_h(n) + Tb_ad*rsh*exptauh2  ! top-of-veg Tb_h (incl atm. contrib.)
             Tb_v(n) = Tb_v(n) + Tb_ad*rsv*exptauv2  ! top-of-veg Tb_v (incl atm. contrib.)

             Tb_h(n) = Tb_h(n) * exptau_atm + Tb_au  ! top-of-atmosphere Tb_h
             Tb_v(n) = Tb_v(n) * exptau_atm + Tb_au  ! top-of-atmosphere Tb_v


          case(2,4)

             ! do nothing

          case default

             err_msg = 'unknown RTM_ID (during Pellarin atm corr)'
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

          end select

       else ! snow present, soil frozen, or mwRTM params not available

          Tb_h(n) = nodata_generic
          Tb_v(n) = nodata_generic

       endif

    end do

    !if (logit) write(logunit,*) 'exiting mwRTM_get_Tb.'

  end subroutine mwRTM_get_Tb


  ! **********************************************************************

  subroutine atmpellarin( Z, tair, costheta, tau_atm, tb_ad, tb_au)

    ! opacity and brightness temperature of atmosphere for low-freq microwave
    ! (up to and incl. L-band)
    !
    ! GDL, 23nov10
    ! Code based on CMEMv3.0
    !
    ! reichle, 16 May 2011: included in LDASsa
    !
    ! ---------------------------------------------------------

    implicit none

    ! arguments

    real, intent(in)  :: costheta    ! cosine of incidence angle      [dim-less]
    real, intent(in)  :: Z           ! elevation above sea-level      [km]
    real, intent(in)  :: tair        ! surface air temperature        [K]

    real, intent(out) :: tau_atm     ! atmospheric opacity            [dim-less]
    real, intent(out) :: tb_ad       ! downwelling atm Tb             [K]
    real, intent(out) :: tb_au       ! upwelling   atm Tb             [K]

    ! local variables

    real              :: GOSSAT, TAEQ

    !---------------------------------------------------------------------
    !
    ! 1. Zenith atmospheric opacity

    tau_atm = exp( -3.926 - 0.2211 * Z - 0.00369 *tair)


    ! 2. Calculate up- and downward atmospheric radiation

    GOSSAT  = exp(-tau_atm/ costheta)

    TAEQ    = exp( 4.927 + 0.002195 * tair)

    tb_ad   = TAEQ*(1. - GOSSAT) + Tb_sky * GOSSAT

    tb_au   = TAEQ*(1. - GOSSAT)

  end subroutine atmpellarin


  ! ************************************************************

   subroutine dielwang( FREQ, WC, TS, wt, wp, poros, sand, clay, eps)

     ! GDL, 28Mar11
     ! Code adapted from CMEMv3.0
     !
     ! reichle, 16 May 2011: included in LDASsa
     !
     !
     ! Purpose :
     !   Calculate the dielectric constant of a wet soil
     !   Developed and validated for 1.4 and 5 GHz.
     !
     ! Reference:
     !  Wang and Schmugge, 1980: An empirical model for the
     !    complex dielectric permittivity of soils as a function of water
     !    content. IEEE Trans. Geosci. Rem. Sens., GE-18, No. 4, 288-295.
     !
     !---------------------------------------------------------------------------

     implicit none

     real, intent(in)     :: FREQ  ! microwave frequency              [Hz]
     real, intent(in)     :: TS    ! soil temperature                 [deg C]

     real, intent(in)     :: WC    ! volumetric soil water content    [m3/m3]

     real, intent(in)     :: wt    ! transition soil moisture         [m3/m3]
     real, intent(in)     :: wp    ! wilting point                    [m3/m3]
     real, intent(in)     :: poros ! porosity
     real, intent(in)     :: sand  ! sand fraction                    [0-1]
     real, intent(in)     :: clay  ! clay fraction                    [0-1]

     complex, intent(out) :: eps   ! dielectric constant of soil-water mixture [dim-less]

     ! -------------------------------------
     !
     ! local variables

     complex, parameter :: j = (0. , 1. )

     real     :: gamma      ! fitting parameter
     real     :: ecl        ! conductivity loss
     complex  :: ew         ! dielectric constant of water
     complex  :: ex         ! dielectric constant of the initially absorbed water

     real     :: alpha

     !---------------------------------------------------------------------------
     !
     ! 0. Compute dielectric constant of free water
     !
     !    assume soil salinity = 0

     call DIEL_WAT( 2, 2, TS, 0., FREQ, clay, sand, poros, wc, ew)

     !---------------------------------------------------------------------------
     !
     ! 1. Calculate dielectric constant of soil-water mixture

     gamma  = -0.57 * wp + 0.481

     ! wt = 0.49 * wp + 0.165     transition SM parameter from Wang and Schmugge 1980;
     !                            note typo in De Lannoy et al 2013 (0.48 instead of 0.49)


     IF (wc <= wt) THEN

        ex  = diel_ice + (ew-diel_ice)*(wc/wt)*gamma

        eps = wc*ex + (poros-wc)*diel_air + (1.-poros)*diel_rock

     ELSE

        ex  = diel_ice + (ew-diel_ice)*gamma

        eps = wt*ex + (wc-wt)*ew + (poros-wc)*diel_air + (1.-poros)*diel_rock

     ENDIF

     !---------------------------------------------------------------------------
     !
     ! 2. add conductivity loss (Wang dielectric model)

     if (FREQ > 2.5e9) then

        alpha = 0.

     else

        alpha = min( 100.*wp, 26.)

     end if

     ecl = alpha * wc**2.

     eps = eps + j * ecl

   end subroutine dielwang


  ! ************************************************************

   SUBROUTINE DIEL_WAT( medium, isal, T, sal, freq, clay, sand, poros, wc, ew)

     ! adapted from CMEM3.0, GDL - 23May11
     !
     ! included in LDASsa, reichle - 2 Jun 2011

     ! Purpose :
     !   Calculate dielectric constant of water in three different media :
     !   pure water, sea water, soil water

     ! Reference:
     !  Dielectric constant of pure water
     !   Ulaby p 2020
     !  Dielectric constant of saline water
     !   1) Stogryn, A. (1971): Equations for calculating the dielectric constant of
     !    saline water, IEEE Transactions on Microwave Theory and Techniques,
     !    Vol. MTT-19, 733-736.
     !   2) Klein, L. A. and C. T. Swift (1977): An improved model
     !    for the dielectric constant of sea water at microwave
     !    frequencies, IEEE Transactions on  Antennas and Propagation,
     !    Vol. AP-25, No. 1, 104-111.
     !  Dielectric constant of soil water
     !   1) Dobson '85. Modified Debye expression
     !         Stern_Gouy double layer theory
     !   2) Ulaby p 2024

     ! Interface :
     !  medium = pure water(0) sea water(1) soil water(2)
     !  isal = Stogryn (1) Klein and Swift (2)

     ! local variables :
     !  N : normality from salinity (Stogryn, modified by Klein and Swift 1977)
     !  T : temperature of water (C)
     !   ew  : dielectric constant of water
     !  sal : water salinity (psu = ppt(weight) )
     !  eps_w0 : static dielectric constant of pure water (Klein and Swift 1977)
     !  eps_sw0 : Static dielectric constant of soil water
     !  tau_w : relaxation time of pure water (stogryn, 1970)
     !  tau_sw : relaxation time of saline water
     !  sigma : ionic conductivity
     !  sigma_eff : effective conductivity of water (S/m)
     !---------------------------------------------------------------------------

     IMPLICIT NONE

     INTEGER, intent(in)  :: medium ! 0=pure water, 1=sea water, 2=soil water
     INTEGER, intent(in)  :: isal   ! 1=Stogryn, 2=Klein and Swift
     REAL,    intent(in)  :: T      ! temperature of water                     [C]
     REAL,    intent(in)  :: sal    ! salinity of water                        [PSU]
     REAL,    intent(in)  :: freq   ! microwave frequency                      [Hz]
     REAL,    intent(in)  :: poros  ! porosity                                 [m3/m3]
     real,    intent(in)  :: sand   ! sand fraction                            [0-1]
     real,    intent(in)  :: clay   ! clay fraction                            [0-1]
     real,    intent(in)  :: wc     ! soil moisture                            [m3/m3]

     COMPLEX, intent(out) :: ew     ! dielec. const. of water at given freq, T, sal, etc.

     ! local variables

     complex, parameter :: j = (0. , 1. )

     real :: rho_b

     REAL :: N, omega, wc_c
     REAL :: sigma_eff
     REAL :: tau_w, tau_sw
     REAL :: eps_w0, eps_sw0, a, bb

     character(len=*), parameter :: Iam = 'DIEL_WAT'
     character(len=400) :: err_msg

     !---------------------------------------------------------------------------

     tau_w = 1.768e-11 + T*(-6.068e-13  + T*(1.104e-14 - T*8.111e-17 ))

     ! tau_w = 1.768e-11  - 6.068e-13 * T  + 1.104e-14 * T**2  - 8.111e-17 * T**3
     !
     ! same as:
     !
     ! tau_w = 1./(2.*pi) * (1.1109e-10 - 3.824e-12 * T + 6.938e-14 * T**2  &
     !                       - 5.096e-16 * T**3)

     omega = 2.0 * MAPL_PI * freq

     rho_b = (1.-poros)*rho_soil         ! soil bulk density [g/cm3]


     SELECT CASE (isal)

     CASE ( 1 ) ! Stogryn (1971)

        N = 0.9141 * sal * (1.707e-2 + sal*(1.205e-5 + sal*4.058e-9))

        eps_sw0 = 87.74 + T*(-0.4008 + T*(9.398e-4 + T*1.410e-6))

        a = 1. + N*(-0.2551  + N*(5.151e-2 - N*6.889e-3))

        eps_sw0 = eps_sw0 * a

        bb = 1. + N*(-0.04896 + 0.1463e-2*T + N*(-0.02967 + N*5.644e-3))

        tau_sw  = tau_w * bb

     CASE ( 2 ) ! Klein and Swift (1977)

        eps_sw0 = 87.134 + T*(-1.949e-1 + T*(-1.276e-2 + T*2.491e-4))

        a = 1. + sal*( 1.613e-5*T -3.656e-3 + sal*(3.210e-5 - sal*4.232e-7))

        eps_sw0 = eps_sw0 * a

        bb = 1. + sal*(2.282e-5*T - 7.638e-4 + sal*(-7.760e-6 + sal*1.105e-8))

        tau_sw  = tau_w * bb

     END SELECT


     SELECT CASE (medium)

     CASE ( 0 ) ! pure water

        eps_w0 = 88.045 + T*(-0.4147 + T*(6.295e-4 + T*1.075e-5))

        ew = diel_watinf + (eps_w0 - diel_watinf) / (1. - j * omega * tau_w)

     CASE ( 1 ) ! sea water

        err_msg = 'medium=1 (sea water) not implemented'
        call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

     CASE ( 2 ) ! soil water

        ! changed units of sand, clay to [0-1], reichle, 2 Jun 2011

        ! Avoid negative sigma_eff for very sandy soils with low bulk densities.

        sigma_eff = max( 0., -1.645 + 1.939*rho_b - 2.256*sand + 1.594*clay)

        ! Modified Debye expression, Dobson '85

        wc_c = MAX(0.001, wc)  ! to avoid dividing by zero

        ew = diel_watinf + (eps_sw0 - diel_watinf) / (1. - j * omega * tau_sw)  &
             + j * sigma_eff / (omega * eps_0) * (rho_soil - rho_b) / (rho_soil * wc_c)

     END SELECT

   END SUBROUTINE DIEL_WAT

   ! **********************************************************************

   SUBROUTINE MIRONOV( freq, mv, clayfrac, er_r )

     ! Soil dielectric mixing model by Mironov et al IEEE TGRS 2009, doi:10.1109/TGRS.2008.2011631
     !
     !  8 May 2023: Implementation taken from SMAP L2_SM_P retrieval model by qliu; clean-up by reichle

     IMPLICIT NONE

     REAL,    INTENT(IN)  :: freq        ! microwave frequency  [Hz]
     REAL,    INTENT(IN)  :: mv          ! volumetric soil water content [m3/m3]
     REAL,    INTENT(IN)  :: clayfrac    ! clay fraction  [0-1]

     COMPLEX, INTENT(OUT) :: er_r        ! complex dielectric constant of moist soil

     ! ------------------------------

     REAL                 :: f, C, mvt
     REAL                 :: nd, nb, nu, nm
     REAL                 :: kd, kb, ku, km
     REAL                 :: eps0b, taub, sigb, epsb_real, epsb_imag
     REAL                 :: eps0u, tauu, sigu, epsu_real, epsu_imag
     REAL                 :: er_r_real, er_r_imag

     REAL                 :: tmp2PIf, tmpeps0

     REAL                 :: tmptaub, tmptaub2plus1, tmpdiffepsb
     REAL                 :: tmptauu, tmptauu2plus1, tmpdiffepsu

     REAL                 :: tmpreal, tmprealb, tmprealu

     ! --------------------------------------------------------------------------------------

     f      = freq                                                        ! Section IV
     C      = clayfrac*100                                                ! Section VI

     !! Mironov's regression expressions based on Curtis, Dobson, and Hallikainen datasets
     !!
     !! mvt    : max bound water fraction
     !!
     !! eps(*) : dielectric  constant (real part) and loss factor (imaginary part)
     !! n(*)   : refractive index
     !! k(*)   : normalized attenuation coefficient
     !! tau(*) : relaxation constant
     !! sig(*) : conductivity
     !!
     !!   m: moist soil
     !!   d: dry soil
     !!   b: bound soil water (BSW)
     !!   u: unbound (free) soil water (FSW)

     nd    = 1.634     -  0.539e-2         * C +  0.2748e-4 * C**2        ! Eqn 17
     kd    = 0.03952   -  0.04038e-2       * C                            ! Eqn 18
     mvt   = 0.02863   +  0.30673e-2       * C                            ! Eqn 19
     eps0b = 79.8      - 85.4e-2           * C + 32.7e-4    * C**2        ! Eqn 20
     taub  = 1.062e-11 +  3.450e-12 * 1e-2 * C                            ! Eqn 21
     sigb  = 0.3112    +  0.467e-2         * C                            ! Eqn 22
     sigu  = 0.3631    +  1.217e-2         * C                            ! Eqn 23
     eps0u = 100.                                                         ! Eqn 24
     tauu  = 8.5e-12                                                      ! Eqn 25

     !! Debye relaxation equations for water as a function of frequency   ! Eqn 16

     tmp2PIf       = 2.*MAPL_PI*f

     tmpeps0       = tmp2PIf*eps_0

     tmptaub       = tmp2PIf*taub
     tmptauu       = tmp2PIf*tauu

     tmptaub2plus1 = 1. + tmptaub**2
     tmptauu2plus1 = 1. + tmptauu**2

     tmpdiffepsb   = eps0b - diel_watinf
     tmpdiffepsu   = eps0u - diel_watinf

     epsb_real = diel_watinf + tmpdiffepsb/tmptaub2plus1
     epsb_imag =               tmpdiffepsb/tmptaub2plus1 * tmptaub + sigb/tmpeps0
     epsu_real = diel_watinf + tmpdiffepsu/tmptauu2plus1
     epsu_imag =               tmpdiffepsu/tmptauu2plus1 * tmptauu + sigu/tmpeps0

     !! Refractive indices and normalized attenuation coefficients

     tmpreal  = 1/sqrt(2.0)

     tmprealb = sqrt( epsb_real**2 + epsb_imag**2 )
     tmprealu = sqrt( epsu_real**2 + epsu_imag**2 )

     nb = tmpreal * sqrt( tmprealb + epsb_real )                          ! Eqn 14
     kb = tmpreal * sqrt( tmprealb - epsb_real )                          ! Eqn 15
     nu = tmpreal * sqrt( tmprealu + epsu_real )                          ! Eqn 14
     ku = tmpreal * sqrt( tmprealu - epsu_real )                          ! Eqn 15

     IF (mv <= mvt) THEN

        nm        = nd + (nb - 1.) * mv                                   ! Eqn 12
        km        = kd +  kb       * mv                                   ! Eqn 13

     ELSE

        nm        = nd + (nb - 1.) * mvt + (nu - 1.) * (mv - mvt)         ! Eqn 12
        km        = kd +  kb       * mvt +  ku       * (mv - mvt)         ! Eqn 13

     ENDIF

     ! complex dielectric constant of moist soil

     er_r_real = nm**2 - km**2                                            ! Eqn 11
     er_r_imag = 2. * nm * km                                             ! Eqn 11

     er_r = CMPLX(er_r_real,er_r_imag)

   END SUBROUTINE mironov

   ! **********************************************************************

end module mwRTM_routines

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if 0

! driver routines for testing

program test_mwRTM_types

  use mwRTM_routines

  implicit none

  type(mwRTM_param_type) :: mwRTM_param

  mwRTM_param = -9999.

  write (*,*) mwRTM_param

end program test_mwRTM_types

#endif

! ========================== EOF ==================================

