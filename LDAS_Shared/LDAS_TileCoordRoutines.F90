
! this file contains types and subroutines for tile coordinates and domain 
!
! reichle, 26 Jan 2005
! reichle,  14 Apr 2006 - split tile_coord.F90 into 2 files to avoid 
!                         having more than one module per file
! reichle,   5 Apr 2013 - added EASEv2 grid, minimal change to max lat/lon for EASE (v1)  
! amfox,     6 Nov 2026 - added cube-sphere grid ("c3") case to get_tile_num_from_latlon()
! amfox,     8 Nov 2026 - map lat/lon to cube-sphere grid cells with MAPL_GetGlobalHorzIJIndex()
!                         and the ESMF grid (get_ij_ind_from_latlon_c3(), init_tile_grid_g_ESMF_c3())
!
! ========================================================================

module LDAS_TileCoordRoutines

  use LDAS_TileCoordType,               ONLY:     &
       tile_coord_type,                           &
       grid_def_type,                             &
       init_grid_def_type,                        &
       io_grid_def_type,                          &
       io_tile_coord_type

  use LDAS_ensdrv_Globals,              ONLY:     &
       logit,                                     &
       logunit,                                   &
       nodata_generic,                            &
       nodata_tol_generic
  
  use MAPL_ConstantsMod,                ONLY:     &
       MAPL_RADIUS,                               &   ! Earth radius
       MAPL_PI,                                   &
       MAPL_PI_R8
  
  use MAPL,                             ONLY:     &
       MAPL_ease_convert,                         &
       MAPL_ease_extent,                          &
       MAPL_GridGet,                              &
       MAPL_GetGlobalHorzIJIndex

  use ESMF,                             ONLY:     &
       ESMF_Grid,                                 &
       ESMF_AttributeGet

  use LDAS_ExceptionsMod,               ONLY:     &
       ldas_abort,                                &
       LDAS_GENERIC_ERROR

  use LDAS_ensdrv_functions,            ONLY:     &
       get_io_filename
  
  implicit none

  private
  
  public :: get_minExtent_grid
  public :: get_number_of_tiles_in_cell_ij
  public :: get_tile_num_in_cell_ij
  public :: get_tile_num_in_ellipse
  public :: get_tile_num_from_latlon
  public :: get_ij_ind_from_latlon
  public :: get_ij_ind_from_latlon_c3
  public :: init_tile_grid_g_ESMF_c3
  public :: tile2grid_simple
  public :: grid2tile
  public :: tile_mask_grid
  public :: LDAS_create_grid_g
  public :: io_domain_files

  character(10) :: tmpstring10
  character(40) :: tmpstring40

  interface grid2tile
    module procedure grid2tile_real, grid2tile_real8
  end interface grid2tile
   
contains
  
  ! **********************************************************************  
  
  subroutine io_domain_files( action, work_path, exp_id,               &
       N_cat_domain, d2g, tile_coord, tile_grid_g, tile_grid_d, rc )
    
    ! reichle, 23 July 2010
    ! reichle,  7 Jan  2014 - changed tile_coord and tile_grids I/O to binary
    ! reichle,  3 Aug  2020 - moved here from LDAS_ensdrv_init_routines.F90
    
    implicit none
    
    character,                          intent(in)    :: action
    
    character(*),                       intent(in)    :: work_path
    character(*),                       intent(in)    :: exp_id

    integer,                             intent(inout) :: N_cat_domain
    
    integer,               dimension(:), pointer       :: d2g

    type(tile_coord_type), dimension(:), pointer       :: tile_coord  ! inout

    type(grid_def_type),                 intent(inout) :: tile_grid_g
    type(grid_def_type),                 intent(inout) :: tile_grid_d

    integer,                             intent(  out) :: rc

    ! local 

    integer, parameter  :: unitnumber = 10
    
    integer             :: n, istat
    
    logical             :: writing
    
    character(300)      :: fname
    character( 40)      :: file_tag, dir_name, file_ext, tmp_action, tmp_status

    character(len=*), parameter :: Iam = 'io_domain_files'
    character(len=400) :: err_msg

    ! -----------------------------------------------------------    
    !
    ! read or write? 
    
    select case (action)
       
    case ('w','W')
       
       tmp_action = 'write'
       tmp_status = 'unknown'

       writing    = .true.

    case ('r','R')
       
       tmp_action = 'read'
       tmp_status = 'old'
       
       writing    = .false.
       
    case default
       
       err_msg = 'io_domain_files: unknown action ' // action
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end select
    
    ! -----------------------------------------------------------    
    !
    !  *.ldas_domain*txt file
    
    file_tag = 'ldas_domain'
    dir_name = 'rc_out'
    file_ext = '.txt'
 
    fname = get_io_filename( trim(work_path), trim(exp_id), file_tag, &
         dir_name=dir_name, option=1, file_ext=file_ext )
 
    if (logit) write (logunit,'(400A)') '  ' // trim(tmp_action) // ' ' // trim(fname)

    open(unitnumber, file=fname, form='formatted', action=trim(tmp_action), &
         iostat=istat, status=trim(tmp_status))
    
    if (istat/=0) then
       
       if (logit) write (logunit,*) 'cannot open for ', trim(tmp_action)
       if (logit) write (logunit,*) 

       rc = 1
       
       return
       
    else
       
       if (writing) then
          
          write (unitnumber,*) N_cat_domain
          
          do n=1,N_cat_domain
             write (unitnumber,*) tile_coord(n)%tile_id, d2g(n)
          end do
          
       else
          
          read  (unitnumber,*) N_cat_domain

          allocate(tile_coord(N_cat_domain))
          allocate(d2g(       N_cat_domain))

          do n=1,N_cat_domain
             read  (unitnumber,*) tile_coord(n)%tile_id, d2g(n)
          end do
          
       end if
       
       close (unitnumber,status='keep')
       
       if (logit) write (logunit,*) 'done with ', trim(tmp_action)
       if (logit) write (logunit,*) 
       
    end if
    
    ! -----------------------------------------------------------    
    !
    ! *.ldas_tile_coord.txt file
    
    file_tag = 'ldas_tilecoord'
    dir_name = 'rc_out'
    file_ext = '.bin'
    
    fname = get_io_filename( work_path, exp_id, file_tag, &
         dir_name=dir_name, option=1, file_ext=file_ext )
    
    if (logit) write (logunit,'(400A)') ' ' // trim(tmp_action) // ' ' // trim(fname)

    open(unitnumber, file=fname, form='unformatted', action=trim(tmp_action), &
         iostat=istat, status=trim(tmp_status))
    
    if (istat/=0) then
       
       if (logit) write (logunit,*) 'cannot open for ', trim(tmp_action)
       if (logit) write (logunit,*) 

       rc = 2
       
       return
       
    else

       call io_tile_coord_type( action, unitnumber, N_cat_domain, tile_coord )
       
       close (unitnumber,status='keep')
       
       if (logit) write (logunit,*) 'done with ', trim(tmp_action)
       if (logit) write (logunit,*) 
       
    end if

    ! -----------------------------------------------------------    
    !
    ! *.ldas_tilegrids.txt file
    
    file_tag = 'ldas_tilegrids'
    dir_name = 'rc_out'
    file_ext = '.bin'
    
    fname = get_io_filename( work_path, exp_id, file_tag, &
         dir_name=dir_name, option=1, file_ext=file_ext )

    if (logit) write (logunit,'(400A)') ' ' // trim(tmp_action) // ' ' // trim(fname)
    
    open(unitnumber, file=fname, form='unformatted', action=trim(tmp_action), &
         iostat=istat, status=trim(tmp_status))
    
    if (istat/=0) then
       
       if (logit) write (logunit,*) 'cannot open for ', trim(tmp_action)
       if (logit) write (logunit,*)        

       rc = 3
       
       return
       
    else
       
       ! read/write 'tile_grid_g'
       call io_grid_def_type( action, unitnumber, tile_grid_g )
       
       ! read/write 'tile_grid_d'
       call io_grid_def_type( action, unitnumber, tile_grid_d )
       
       close (unitnumber,status='keep')
       
       if (logit) write (logunit,*) 'done with ', trim(tmp_action)
       if (logit) write (logunit,*) 

    end if
    
    ! ------------------------------------------------------------------------
    
    rc = 0    ! successful read or write of all files
    
  end subroutine io_domain_files

  ! **********************************************************************

  subroutine LDAS_create_grid_g( gridname, n_lon, n_lat,       &
       tile_grid, i_indg_offset, j_indg_offset, cell_area)
    
    ! inputs:
    !  gridname, n_lon, n_lat
    !
    ! inouts:
    !  tile_grid   : parameters of tile definition grid  
    ! 
    ! outputs:
    !  offsets
    !  cell_area  [m^2]   (optional, for EASE  grids only)
    
    implicit none
    
    character(*),                  intent(in)    :: gridname
    integer,                       intent(in)    :: n_lon, n_lat
    type(grid_def_type),           intent(inout) :: tile_grid
    integer,                       intent(out)   :: i_indg_offset, j_indg_offset
    real,                optional, intent(out)   :: cell_area
    
    ! locals

    real    :: ease_cell_area
    logical :: date_line_on_center, pole_on_center
    logical :: ease_grid, c3_grid, latlon_grid
    logical :: file_exist
    integer :: k, rows, cols
    character(len=*),  parameter :: Iam = 'create global ldas_grid '
    character(len=400)           :: err_msg
    
    ! initialize all fields to no-data values
    
    i_indg_offset = 0
    j_indg_offset = 0
    
    call init_grid_def_type(tile_grid)
    
    tile_grid%N_lon = N_lon
    tile_grid%N_lat = N_lat
    
    tile_grid%i_offg = 0  ! tile_grid refers to *global* grid
    tile_grid%j_offg = 0  ! tile_grid refers to *global* grid
    
    date_line_on_center = .false.
    pole_on_center      = .false.
    ease_grid           = .false.
    c3_grid             = .false.
    latlon_grid         = .true.
    
    if (index(gridname,"DC") /=0) then
       date_line_on_center = .true.
    endif
    
    if (index(gridname,"PC") /=0) then
       pole_on_center      = .true.
    endif
    
    if( index(gridname,"FV") /=0 ) then
       pole_on_center      = .true.
    endif
    
    if (index(gridname,"EASEv") /=0) then
       ease_grid      = .true.
       latlon_grid    = .false.
    endif
    
    if (index(gridname,"CF") /=0) then
       c3_grid        = .true.
       latlon_grid    = .false.
    endif
    
    ! special cases , inconsistent of naming
    ! find out whether date line is on edge or through center of grid cell
    if( index(gridname,"FV_380x180") /=0) then
       pole_on_center      = .false.
    endif
    
    ! Weiyuan Note, we should fix the tile file and the naming 
    if( index(gridname,"PE_720x360_DE") /=0) then
       i_indg_offset = 110
       j_indg_offset = 230
    endif
    if( index(gridname,"PE_2880x1440_DE") /=0) then
       i_indg_offset = 440
       j_indg_offset = 920
    endif
    
    ! ----------------
    
    if (ease_grid) then
       
       ! gridname may be EASEv2-M36 or EASEv2_M36 (to be cleaned up later)

       k = index(gridname, 'EASEv')
       
       if (k == 0) then   
          err_msg = 'unknown EASE grid tile defs, gridname = ' // trim( gridname)
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       end if
       
       tile_grid%gridtype = trim(gridname(k:))
       
       tile_grid%ind_base = 0
       
       ! global cylindrical EASE grid 
       
       tile_grid%i_dir  = +1
       tile_grid%j_dir  = -1
       
       call MAPL_ease_extent (              &
            gridname, cols, rows,           &
            cell_area = ease_cell_area,     &             ! [m^2]
            ll_lon    = tile_grid%ll_lon,   &
            ll_lat    = tile_grid%ll_lat,   &
            ur_lon    = tile_grid%ur_lon,   &
            ur_lat    = tile_grid%ur_lat    &
            )

       tile_grid%dlon   = 360./real(tile_grid%N_lon) 
       tile_grid%dlat   = (tile_grid%ur_lat-tile_grid%ll_lat)/real(tile_grid%N_lat) ! *avg* dlat!

       if(present(cell_area)) then
          cell_area=ease_cell_area
       endif
       
    endif ! EASE grid
    
    if (latlon_grid) then !  regular LatLon grid
       
       tile_grid%gridtype = 'LatLon'
       
       tile_grid%ind_base = 1
       
       tile_grid%dlon = 360./real(tile_grid%N_lon)
       
       tile_grid%i_dir = +1
       tile_grid%j_dir = +1
       
       if (pole_on_center) then
          
          tile_grid%dlat   = 180./real(tile_grid%N_lat-1)
          
          tile_grid%ll_lat = -90. - tile_grid%dlat/2.
          tile_grid%ur_lat =  90. + tile_grid%dlat/2.
          
       else
          
          tile_grid%dlat = 180./real(tile_grid%N_lat)
          
          tile_grid%ll_lat = -90. 
          tile_grid%ur_lat =  90. 
          
       end if
       
       if (date_line_on_center) then
          
          tile_grid%ll_lon = -180. - tile_grid%dlon/2.
          tile_grid%ur_lon =  180. - tile_grid%dlon/2.  ! fixed 20 sep 2010, reichle
          
       else
          
          tile_grid%ll_lon = -180. 
          tile_grid%ur_lon =  180. 
          
       end if
       
    end if ! lat lon grid
    
    if( c3_grid) then
       
       tile_grid%gridtype = 'c3'
       tile_grid%ind_base =    1
       tile_grid%i_dir    =   +1
       tile_grid%j_dir    =   +1
       tile_grid%ll_lon   = -180. 
       tile_grid%ur_lon   =  180. 
       tile_grid%ll_lat   =  -90. 
       tile_grid%ur_lat   =   90.
       
       ! dlon and dlat are approximate!
       tile_grid%dlon     = 360./real(4*tile_grid%N_lon)
       tile_grid%dlat     = tile_grid%dlon
       
    endif
    
  end subroutine LDAS_create_grid_g
  
  ! *******************************************************************
  
  function get_minExtent_grid( N_tile, tc_i_indg, tc_j_indg, tc_minlon, tc_minlat, tc_maxlon, tc_maxlat, &
       tile_grid_g) result( tile_grid )     
    
    ! get tile_grid with smallest extent for given set of tiles in tile_coord on (global) tile_grid_g
    !
    ! make sure to pass in consistent tile_coord (tc) and tile_grid_g inputs:
    !
    !   iff tile_grid_g is the grid that is associated with the tile space definition,
    !   then input tc_[x]_indg must be tile_coord%[x]_indg
    !
    !   iff tile_grid_g is the grid that supports efficient mapping for perts and the EnKF analysis,
    !   then input tc_[x]_indg must be tile_coord%pert_[x]_indg
    !
    ! reichle, 20 June 2012 -- moved from within domain_setup() 
    !                           for re-use in get_obs_pred() 
    !
    ! reichle+wjiang, 19 Aug 2020 -- changed interface to generically accommodate use of
    !                                tile_coord%[x]_indg or tile_coord%pert_[x]_indg
    !
    ! -------------------------------------------------------------------
    
    integer,               intent(in)                    :: N_tile

    integer,               intent(in), dimension(N_tile) :: tc_i_indg, tc_j_indg
    real,                  intent(in), dimension(N_tile) :: tc_minlon, tc_minlat
    real,                  intent(in), dimension(N_tile) :: tc_maxlon, tc_maxlat
    
    type(grid_def_type),   intent(in)                    :: tile_grid_g
    type(grid_def_type)                                  :: tile_grid
    
    ! local variables
    
    integer :: n
    
    real    :: this_minlon, this_minlat, this_maxlon, this_maxlat
    
    real    :: min_min_lon, min_min_lat, max_max_lon, max_max_lat
    
    integer :: ind_i_min, ind_i_max, ind_j_min, ind_j_max
    
    integer :: this_i_indg, this_j_indg
    integer , allocatable :: i_indg_(:), j_indg_(:)

    character(len=*), parameter :: Iam = 'get_minExtent_grid'
    character(len=400) :: err_msg

    ! -------------------------------------------------
    
    min_min_lon =  180.                  ! initialize
    min_min_lat =   90.
    max_max_lon = -180.
    max_max_lat =  -90.
    
    ind_i_min = tile_grid_g%N_lon+1      ! initialize
    ind_j_min = tile_grid_g%N_lat+1
    ind_i_max = -1
    ind_j_max = -1

    ! THIS COMMENT SEEMS OUTDATED (reichle, 2 Aug 2020)
    ! for c3 grid, only get the ll_,ur_ lat and lon, the index is meaningless;
    ! it will be used in creating the lat_lon pert_grid
    
    if(index(tile_grid_g%gridtype,"c3") /=0) then

       ! for cube-sphere grids, do NOT zoom in 
       
       tile_grid=tile_grid_g
       
       return
       
    endif   

    if (N_tile == 0) then
      tile_grid=tile_grid_g
      tile_grid%n_lon = 0
      tile_grid%n_lat = 0
      return
    endif

    do n=1,N_tile
       
       this_minlon  = tc_minlon(n)
       this_minlat  = tc_minlat(n)
       this_maxlon  = tc_maxlon(n)
       this_maxlat  = tc_maxlat(n)
       
       this_i_indg  = tc_i_indg(n)
       this_j_indg  = tc_j_indg(n)
       
       min_min_lon = min( min_min_lon, this_minlon)
       min_min_lat = min( min_min_lat, this_minlat)
       max_max_lon = max( max_max_lon, this_maxlon)
       max_max_lat = max( max_max_lat, this_maxlat)
       
       ind_i_min   = min( ind_i_min,   this_i_indg)
       ind_j_min   = min( ind_j_min,   this_j_indg)
       ind_i_max   = max( ind_i_max,   this_i_indg)
       ind_j_max   = max( ind_j_max,   this_j_indg)
       
    end do
    
    ! assemble tile_grid (revised 20 Sep 2010, reichle)
       
    tile_grid%N_lon  = ind_i_max - ind_i_min + 1
    tile_grid%N_lat  = ind_j_max - ind_j_min + 1
    
    tile_grid%i_offg = ind_i_min - tile_grid_g%ind_base
    tile_grid%j_offg = ind_j_min - tile_grid_g%ind_base
       
    tile_grid%dlon   = tile_grid_g%dlon
       
    if     (index(tile_grid_g%gridtype,'LatLon')/=0) then
       
       tile_grid%dlat   = tile_grid_g%dlat
       
       tile_grid%ll_lon = tile_grid_g%ll_lon + real(tile_grid%i_offg)*tile_grid%dlon
       tile_grid%ll_lat = tile_grid_g%ll_lat + real(tile_grid%j_offg)*tile_grid%dlat
       
       tile_grid%ur_lon = tile_grid%ll_lon + real(tile_grid%N_lon)*tile_grid%dlon
       tile_grid%ur_lat = tile_grid%ll_lat + real(tile_grid%N_lat)*tile_grid%dlat
       
    elseif ( index(tile_grid_g%gridtype,'EASEv')  /=0 ) then
       
       ! *average* dlat over the domain
       
       tile_grid%dlat   = (max_max_lat - min_min_lat)/real(tile_grid%N_lat)
       
       tile_grid%ll_lon = min_min_lon
       tile_grid%ll_lat = min_min_lat
       
       tile_grid%ur_lon = max_max_lon
       tile_grid%ur_lat = max_max_lat
       
    else
       
       err_msg = 'not yet implemented for ' // tile_grid_g%gridtype
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
       
    end if
    
    tile_grid%gridtype = tile_grid_g%gridtype
    
    tile_grid%ind_base = tile_grid_g%ind_base   ! should this be inherited???
    
    tile_grid%i_dir    = tile_grid_g%i_dir
    tile_grid%j_dir    = tile_grid_g%j_dir
    
  end function get_minExtent_grid

  ! **********************************************************************

  subroutine get_number_of_tiles_in_cell_ij( N_tile, tc_i_indg, tc_j_indg, tile_grid, &
       N_tile_in_cell_ij)
    
    ! find out how many tiles are in a given tile definition grid cell
    ! reichle, 22 Jul 2005
    !
    ! make sure to pass in consistent tile_coord (tc) and tile_grid_g inputs:
    !
    !   iff tile_grid_g is the grid that is associated with the tile space definition,
    !   then tc_[x]_indg must be tile_coord%[x]_indg
    !
    !   iff tile_grid_g is the grid that supports efficient mapping for perts and the EnKF analysis,
    !   then tc_[x]_indg must be tile_coord%pert_[x]_indg
    !
    ! wjiang(?) -- split off from LDASsa legacy subroutine get_tile_num_in_cell_ij()
    !
    ! reichle+wjiang, 19 Aug 2020 -- changed interface to generically accommodate use of
    !                                tile_coord%[x]_indg or tile_coord%pert_[x]_indg
    !
    ! ----------------------------------------------------------
    
    implicit none
    
    integer,             intent(in)                    :: N_tile
    
    integer,             intent(in), dimension(N_tile) :: tc_i_indg, tc_j_indg
    
    type(grid_def_type), intent(in)                    :: tile_grid
        
    integer, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(inout) :: &
         N_tile_in_cell_ij
    
    ! locals 
    
    integer :: i, j, k, n, off_i, off_j
    

    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    ! (re-)initialize
    
    N_tile_in_cell_ij = 0
    
    do n=1,N_tile
       
       i = tc_i_indg(n) - off_i
       j = tc_j_indg(n) - off_j
       
       N_tile_in_cell_ij(i,j) = N_tile_in_cell_ij(i,j) + 1
       
    enddo

    if (logit) write (logunit,*) &
        'Maximum number of tiles in tile def grid cell = ', maxval(N_tile_in_cell_ij)
    if (logit) write (logunit,*)

  end subroutine get_number_of_tiles_in_cell_ij

  ! **********************************************************************

  subroutine get_tile_num_in_cell_ij( N_tile, tc_i_indg, tc_j_indg, tile_grid, &
       max_N_tile_in_cells, tile_num_in_cell_ij )
    
    ! find out tile_num in given cells
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domain of interest.  
    ! With these offsets tile2grid() can be used to map from a
    ! subgrid of "tile_grid_g" to tile space
    !
    ! reichle, 22 Jul 2005
    !
    ! make sure to pass in consistent tile_coord (tc) and tile_grid_g inputs:
    !
    !   iff tile_grid_g is the grid that is associated with the tile space definition,
    !   then tc_[x]_indg must be tile_coord%[x]_indg
    !
    !   iff tile_grid_g is the grid that supports efficient mapping for perts and the EnKF analysis,
    !   then tc_[x]_indg must be tile_coord%pert_[x]_indg
    !
    ! wjiang(?) -- split off from LDASsa legacy subroutine get_tile_num_in_cell_ij()
    !
    ! reichle+wjiang, 19 Aug 2020 -- changed interface to generically accommodate use of
    !                                tile_coord%[x]_indg or tile_coord%pert_[x]_indg
    !
    ! ----------------------------------------------------------
    
    implicit none
    
    integer,             intent(in)                    :: N_tile
    
    integer,             intent(in), dimension(N_tile) :: tc_i_indg, tc_j_indg
    
    type(grid_def_type), intent(in)                    :: tile_grid
        
    integer,             intent(in)                    ::  max_N_tile_in_cells
    
    ! the pointer is an output arguments that is allocated here
    
    integer, dimension(:,:,:), pointer :: tile_num_in_cell_ij

    integer, dimension(tile_grid%N_lon,tile_grid%N_lat) :: &
         N_tile_in_cell_ij

    ! locals 
    
    integer, parameter :: nodata = -9999.
    
    integer :: i, j, k, n, off_i, off_j
    
    ! -----------------------------------------------------------------
    !
    ! allocate and initialize pointers if present
    
       
    allocate(tile_num_in_cell_ij(tile_grid%N_lon,tile_grid%N_lat,    &
          max_N_tile_in_cells))
       
    tile_num_in_cell_ij = nodata
       
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    ! (re-)initialize
    
    N_tile_in_cell_ij = 0
    
    do n=1,N_tile
       
       i = tc_i_indg(n) - off_i
       j = tc_j_indg(n) - off_j
       
       N_tile_in_cell_ij(i,j) = N_tile_in_cell_ij(i,j) + 1
       
       k = N_tile_in_cell_ij(i,j)
          
       tile_num_in_cell_ij(i,j,k) = n 
          
    end do
    
  end subroutine get_tile_num_in_cell_ij
  
  ! *******************************************************************
  
  subroutine get_tile_num_from_latlon(N_catd, tile_coord,                      &
       tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij, N_latlon, lat, lon,  &
       tile_num, max_dist_x, max_dist_y, tile_grid_g_ESMF )

    ! bug fix re. "check that lat/lon is inside tile_grid"
    ! - reichle, 2005/11/17
    !
    ! added optional input "max_dist" that permits returning a valid tile_num
    ! even if lat/lon is not within the bounding box of the chosen tile
    ! - reichle, 2008/03/28
    !
    ! distance units are in [deg]; added functionality for latitude-dependent
    !  meridional max distance (max_dist_x)
    ! - reichle, 2014/12/23
    !
    ! added cube-sphere ("c3") tile_grid, see get_tile_num_from_latlon_c3();
    !  requires optional input "tile_grid_g_ESMF" (ESMF grid of cube-sphere tile space)
    ! - amfox, 8 Nov 2026
    
    implicit none
    
    integer,                                             intent(in)  :: N_catd, N_latlon
    
    type(tile_coord_type), dimension(:),                 pointer     :: tile_coord ! input
    
    type(grid_def_type),                                 intent(in)  :: tile_grid    
    
    integer, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in)  :: N_tile_in_cell_ij
    
    integer, dimension(:,:,:),                           pointer     :: tile_num_in_cell_ij ! input    
    
    real,    dimension(N_latlon),                        intent(in)  :: lat, lon
    
    integer, dimension(N_latlon),                        intent(out) :: tile_num
    
    real,    dimension(N_latlon),              optional, intent(in)  :: max_dist_x  ! vector [deg]
    real,                                      optional, intent(in)  :: max_dist_y  ! scalar [deg]
    
    type(ESMF_Grid),                           optional, intent(in)  :: tile_grid_g_ESMF     ! for "c3" tile_grid

    ! local variables
    
    integer, parameter           :: nodata_tilenum = -9999
    
    integer                      :: n, k, i_ind, j_ind, this_tile_num

    logical                      :: outside_bbox, too_far_away
    
    real                         :: tmp_dist,   min_dist
    real                         :: tmp_dist_y, min_dist_y, max_dist_tmp_y
    real                         :: tmp_dist_x, min_dist_x
    real, dimension(N_latlon)    ::                         max_dist_tmp_x 
    
    character(len=*),  parameter :: Iam = 'get_tile_num_from_latlon'
    character(len=400)           :: err_msg

    ! -----------------------------------------------------------
    
    ! initialize
    
    max_dist_tmp_x = 0.
    max_dist_tmp_y = 0.

    if (present(max_dist_x))  max_dist_tmp_x = max_dist_x
    if (present(max_dist_y))  max_dist_tmp_y = max_dist_y
    
    tile_num = nodata_tilenum          ! initialize to negative value
    
    ! cube-sphere tile_grid: tiles do not nest in grid cells of a lat/lon grid,
    !  use bounding boxes of tiles (amfox, 8 Nov 2026)

    if (index(tile_grid%gridtype, 'c3')/=0) then

       if (.not. present(tile_grid_g_ESMF)) then

          err_msg = 'cube-sphere ("c3") tile_grid requires optional input tile_grid_g_ESMF (ESMF grid)'
          call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

       end if

       call get_tile_num_from_latlon_c3( tile_grid_g_ESMF, tile_coord,              &
            tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij, N_latlon, lat, lon,  &
            max_dist_tmp_x, max_dist_tmp_y, tile_num )

       return

    end if

    ! loop through observations
    
    do n=1,N_latlon
       
       ! Make sure lat/lon is *inside* tile_grid, otherwise do nothing.
       ! Do *not* allow obs outside tile_grid because obs perturbations 
       ! routines are not set up to generate perturbations for such obs, 
       ! and obs_pred would be questionable anyway.
       ! - reichle+csdraper, 29 Jan 2016

       if ( tile_grid%ll_lat   < (lat(n)          )  .and.  &
            tile_grid%ll_lon   < (lon(n)          )  .and.  &
            (lat(n)          ) < tile_grid%ur_lat    .and.  &
            (lon(n)          ) < tile_grid%ur_lon              ) then  

          ! min_dist = distance betw lat/lon in question and center-of-mass of
          !            matching tile 
          
          min_dist   = 1.e10    ! initialize   (bug fix, csdraper+reichle, 30 Jun 2015)
          min_dist_x = 1.e10    ! initialize 
          min_dist_y = 1.e10    ! initialize 
          
          ! determine grid cell that contains lat/lon 
          
          call get_ij_ind_from_latlon( tile_grid, lat(n), lon(n), i_ind, j_ind )
          
          ! make sure that i/j_ind is still within bounds 
          ! (works in conjunction with if statement above re. ll/ur_lat/lon)
          
          i_ind = min( max(i_ind, 1), tile_grid%N_lon )
          j_ind = min( max(j_ind, 1), tile_grid%N_lat )
          
          ! map from i_ind, j_ind to tile_num
          
          if   ( index(tile_grid%gridtype, 'EASEv')  /=0 ) then
             
             ! ASSUMPTION: tiles match EASE grid cells exactly
             !             (unless "outside" the domain, eg. water surface)
             
             if     (N_tile_in_cell_ij(i_ind,j_ind)==1) then
                
                tile_num(n)=tile_num_in_cell_ij(i_ind,j_ind,1)

                min_dist_x = abs(lon(n) - tile_coord(tile_num(n))%com_lon)
                min_dist_y = abs(lat(n) - tile_coord(tile_num(n))%com_lat) 
                
             elseif (N_tile_in_cell_ij(i_ind,j_ind)==0) then
                
                ! Do nothing.  If given EASE grid cell is not land, 
                ! tile_num will not change from its initialized value.
                
             else
                
                err_msg = 'something wrong for EASE grid'
                call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
                
             end if
             
          elseif (index(tile_grid%gridtype, 'LatLon')/=0) then
             
             ! Loop through all tiles within grid cell that contain lat/lon and
             ! find minimum distance (Minkowski norm).
             ! If there are no land tiles in given tile definition grid cell, 
             ! tile_num will not change from its initialized value.
             
             do k=1,N_tile_in_cell_ij(i_ind,j_ind) 
                
                this_tile_num = tile_num_in_cell_ij(i_ind,j_ind,k)

                tmp_dist_x = abs(lon(n) - tile_coord(this_tile_num)%com_lon)
                tmp_dist_y = abs(lat(n) - tile_coord(this_tile_num)%com_lat) 

                tmp_dist   = tmp_dist_x + tmp_dist_y    ! Minkowski norm
                
                if (tmp_dist<min_dist) then
                   
                   min_dist    = tmp_dist
                   min_dist_x  = tmp_dist_x
                   min_dist_y  = tmp_dist_y
                   
                   tile_num(n) = this_tile_num
                   
                end if
                
             end do

          else
             
             call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown grid type')
             
          end if
                    

          ! check whether lat/lon is inside bounding box of given tile
          ! (if domain consists of a single tile and tile_grid is coarse,
          !  it is likely that lat/lon is outside bounding box of tile
          !  even if the lat/lon is inside the tile definition grid cell that
          !  contains the given tile) - reichle, 18 Aug 2005
          
          ! return a valid tile_num even if lat/lon is not within the 
          ! bounding box of the chosen tile, as long as lat/lon is no further
          ! from the com_lat/com_lon of the tile than max_dist
          ! - reichle, 2008/03/28
          
          ! this might be important for assimilation of 36km EASE grid obs 
          ! into the 9 km EASE grid tile space (because the lat/lon of
          ! the 36 km EASE grid obs is exactly on the corners of four
          ! 9 km EASE grid cells)
          ! - reichle, 2013/06/21          

          if (tile_num(n)>0) then

             outside_bbox  = (                                          &
                  lon(n) < tile_coord(tile_num(n))%min_lon    .or.      &
                  lon(n) > tile_coord(tile_num(n))%max_lon    .or.      &
                  lat(n) < tile_coord(tile_num(n))%min_lat    .or.      &
                  lat(n) > tile_coord(tile_num(n))%max_lat           )      
             
             too_far_away = (                                           &
                  min_dist_x > max_dist_tmp_x(n)              .or.      &
                  min_dist_y > max_dist_tmp_y                        )
             
             ! keep tile_num unless obs is outside the bounding box *and* too far away
             
             if (outside_bbox .and. too_far_away)  tile_num(n) = nodata_tilenum
             
          end if
          
       end if
       
    end do
        
  end subroutine get_tile_num_from_latlon
  
  ! *******************************************************************
  
  subroutine get_tile_num_from_latlon_c3( tile_grid_g_ESMF, tile_coord,        &
       tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij, N_latlon, lat, lon,  &
       max_dist_x, max_dist_y, tile_num )

    ! find tile_num for lat/lon in cube-sphere ("c3") tile space, see get_tile_num_from_latlon()
    !
    ! tiles do not nest in the cells of a lat/lon grid, so for cube-sphere tile space,
    !  tile_grid must be the cube-sphere grid that defines the tile space (and
    !  N_tile_in_cell_ij, tile_num_in_cell_ij must be obtained from tile_coord%i_indg and
    !  tile_coord%j_indg).  tile_grid_g_ESMF is the ESMF grid of the cube-sphere tile space.
    !
    ! Each tile is inside the cube-sphere grid cell given by tile_coord%i_indg and
    ! tile_coord%j_indg, but a tile can be much smaller than the grid cell, and lat/lon is
    ! often not inside the tile with the nearest center-of-mass.  Therefore, use the
    ! bounding box of the tiles:
    !
    ! 1.) Loop through all tiles within grid cell that contains lat/lon.
    !     Among the tiles whose bounding box contains lat/lon, find the minimum
    !     normalized square distance (dx/hx)**2 + (dy/hy)**2, where hx and hy
    !     are half the width and height of the bounding box.
    !
    ! 2.) If there is no such tile, repeat 1.) for all tiles within the grid cells
    !     that contain (i) the 8 lat/lon points at +/- half a grid cell (approx.) around
    !     lat/lon, because the bounding box of a tile can extend into a
    !     neighboring grid cell, and (ii) the points of a lattice that covers
    !     lat/lon +/- max_dist, with a spacing of at most a quarter grid cell (approx.).
    !     The lattice finds nearly all grid cells within max_dist, but it can miss a grid
    !     cell that overlaps only a small corner of the lat/lon +/- max_dist rectangle
    !     (3 of 83,645 Step 3.) cases in an offline test for CF0360 with FOV 20 and 50 km;
    !     another tile within max_dist was then chosen).  Grid cells that were already
    !     searched are skipped.
    !
    ! 3.) If there is still no such tile, loop through all tiles within the grid
    !     cells identified in 1.) and 2.) and find minimum distance (Minkowski norm)
    !     among the tiles whose center-of-mass is within max_dist_x and max_dist_y
    !     of lat/lon.
    !
    ! The tiles identified in this way are still subject to the bounding box and
    !   max_dist checks of get_tile_num_from_latlon(), which they pass by construction.
    !
    ! If there are no land tiles in these grid cells, tile_num will not change
    ! from its initialized value.
    !
    ! The grid cells are obtained from MAPL_GetGlobalHorzIJIndex() (see
    !  get_ij_ind_from_latlon_c3()), which has a fixed cost per call.  lat/lon are
    !  therefore processed in chunks, with one call for 1.) and one call per batch of at
    !  most N_pts_max points for 2.).  The fixed cost grows with the grid resolution (~4x
    !  per doubling; ~0.1 s per call for C2880), whereas the cost per point is small
    !  (~0.05 us), hence the large chunks and batches.
    !
    ! amfox, 6 Nov 2026
    ! amfox, 8 Nov 2026 - moved here from get_tile_num_from_latlon(); grid cells from
    !                     MAPL_GetGlobalHorzIJIndex() (no copied cube-sphere geometry)
    ! amfox, 8 Nov 2026 - Step 3.) considers only tiles within max_dist (otherwise the
    !                     closest tile could fail the max_dist check while another tile
    !                     within max_dist would pass it)
    ! amfox, 8 Nov 2026 - Step 2.) also searches the grid cells within max_dist (the 8 points
    !                     alone can miss grid cells with tiles within max_dist)
    ! amfox, 9 Nov 2026 - larger N_chunk and N_pts_max (fewer MAPL calls for high-res grids)

    implicit none

    type(ESMF_Grid),                                     intent(in)  :: tile_grid_g_ESMF

    type(tile_coord_type), dimension(:),                 pointer     :: tile_coord ! input

    type(grid_def_type),                                 intent(in)  :: tile_grid

    integer, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in)  :: N_tile_in_cell_ij

    integer, dimension(:,:,:),                           pointer     :: tile_num_in_cell_ij ! input

    integer,                                             intent(in)  :: N_latlon

    real,    dimension(N_latlon),                        intent(in)  :: lat, lon

    real,    dimension(N_latlon),                        intent(in)  :: max_dist_x  ! vector [deg]
    real,                                                intent(in)  :: max_dist_y  ! scalar [deg]

    integer, dimension(N_latlon),                        intent(out) :: tile_num

    ! local variables

    integer, parameter           :: nodata_tilenum = -9999

    integer, parameter           :: N_chunk        = 100000  ! max number of lat/lon per chunk

    integer, parameter           :: N_lattice_max  = 20      ! max lattice intervals on either side
                                                             !  of lat/lon (Step 2.)

    integer, parameter           :: N_cell_max     = 9 + (2*N_lattice_max+1)**2

    integer, parameter           :: N_pts_max      = 1000000 ! max number of points per batch (Step 2.)

    integer                      :: n, nn, k, kk, ii, jj, mm, this_tile_num
    integer                      :: ind_start, ind_end, N_tmp, N_nbr, N_cell
    integer                      :: nn_start, nn_end, N_x, N_y

    logical                      :: outside_bbox, too_far_away

    real                         :: tmp_dist,   min_dist
    real                         :: tmp_dist_x, tmp_dist_y, tmp_dlon, tmp_dlat, tmp_slon, tmp_slat

    integer, dimension(N_cell_max) :: i_cell_c3, j_cell_c3

    integer, dimension(:),   allocatable :: i_ind, j_ind, ind_nbr, N_pts_nbr, i_ind_nbr, j_ind_nbr

    real,    dimension(:),   allocatable :: min_dist_x, min_dist_y, tmp_lat, tmp_lon

    ! -----------------------------------------------------------

    tile_num = nodata_tilenum          ! initialize to negative value

    allocate( i_ind(N_chunk), j_ind(N_chunk), ind_nbr(N_chunk), N_pts_nbr(N_chunk) )
    allocate( min_dist_x(N_chunk), min_dist_y(N_chunk) )
    allocate( i_ind_nbr(N_pts_max), j_ind_nbr(N_pts_max), tmp_lat(N_pts_max), tmp_lon(N_pts_max) )

    ! loop through chunks of observations

    do ind_start=1,N_latlon,N_chunk

       ind_end = min( ind_start+N_chunk-1, N_latlon )

       N_tmp   = ind_end - ind_start + 1

       min_dist_x = 1.e10    ! initialize
       min_dist_y = 1.e10    ! initialize

       ! determine grid cells that contain lat/lon

       call get_ij_ind_from_latlon_c3( tile_grid_g_ESMF, tile_grid, N_tmp,                    &
            lat(ind_start:ind_end), lon(ind_start:ind_end), i_ind(1:N_tmp), j_ind(1:N_tmp) )

       ! Step 1.): tiles within grid cell that contains lat/lon

       N_nbr = 0

       do nn=1,N_tmp

          n = ind_start - 1 + nn

          ! Make sure lat/lon is *inside* tile_grid, otherwise do nothing
          ! (see get_tile_num_from_latlon())

          if ( tile_grid%ll_lat   < (lat(n)          )  .and.  &
               tile_grid%ll_lon   < (lon(n)          )  .and.  &
               (lat(n)          ) < tile_grid%ur_lat    .and.  &
               (lon(n)          ) < tile_grid%ur_lon              ) then

             min_dist = 1.e10    ! initialize

             call get_tile_num_in_bbox_c3( lat(n), lon(n), 1, i_ind(nn:nn), j_ind(nn:nn), &
                  tile_coord, tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij,         &
                  tile_num(n), min_dist, min_dist_x(nn), min_dist_y(nn) )

             if (tile_num(n)==nodata_tilenum) then

                ! no tile found: lat/lon needs Step 2.)

                N_nbr          = N_nbr + 1
                ind_nbr(N_nbr) = nn

             end if

          end if

       end do

       ! Step 2.): add grid cells that contain points around lat/lon:
       !
       !  (i)  8 points at +/- half a grid cell (approx.) around lat/lon
       !       (tile_grid%dlon = 360/(4*N) is the approximate grid spacing, see LDAS_create_grid_g())
       !
       !  (ii) points of a lattice that covers lat/lon +/- max_dist, with a spacing of at most
       !       a quarter grid cell (approx.); the number of lattice intervals on either side of
       !       lat/lon is limited to N_lattice_max (only for very large max_dist or near the poles)
       !
       ! points are processed in batches of at most N_pts_max points (one MAPL call per batch)

       nn_start = 1

       do while (nn_start<=N_nbr)

          ! points for lat/lon ind_nbr(nn_start:nn_end)

          kk     = 0
          nn_end = nn_start - 1

          do nn=nn_start,N_nbr

             n = ind_start - 1 + ind_nbr(nn)

             tmp_dlat = 0.5*tile_grid%dlon
             tmp_dlon = tmp_dlat / max( cos( MAPL_PI/180. * lat(n) ), 0.1 )

             tmp_slat = 0.25*tile_grid%dlon
             tmp_slon = tmp_slat / max( cos( MAPL_PI/180. * lat(n) ), 0.1 )

             N_y = ceiling( min( max_dist_y   /tmp_slat, real(N_lattice_max) ) )
             N_x = ceiling( min( max_dist_x(n)/tmp_slon, real(N_lattice_max) ) )

             N_pts_nbr(nn) = 8 + (2*N_x+1)*(2*N_y+1)

             if (kk+N_pts_nbr(nn) > N_pts_max)  exit     ! batch is full

             ! (i) 8 points at +/- half a grid cell

             do ii=-1,1
                do jj=-1,1

                   if (ii==0 .and. jj==0) cycle

                   kk = kk + 1

                   tmp_lat(kk) = lat(n) + real(jj)*tmp_dlat
                   tmp_lon(kk) = lon(n) + real(ii)*tmp_dlon

                end do
             end do

             ! (ii) lattice over lat/lon +/- max_dist (just lat/lon if max_dist=0)

             do ii=-N_x,N_x
                do jj=-N_y,N_y

                   kk = kk + 1

                   tmp_lat(kk) = lat(n)
                   tmp_lon(kk) = lon(n)

                   if (N_y>0)  tmp_lat(kk) = lat(n) + real(jj)*(max_dist_y   /real(N_y))
                   if (N_x>0)  tmp_lon(kk) = lon(n) + real(ii)*(max_dist_x(n)/real(N_x))

                end do
             end do

             nn_end = nn

          end do

          call get_ij_ind_from_latlon_c3( tile_grid_g_ESMF, tile_grid, kk,               &
               tmp_lat(1:kk), tmp_lon(1:kk), i_ind_nbr(1:kk), j_ind_nbr(1:kk) )

          kk = 0

          do nn=nn_start,nn_end

             n = ind_start - 1 + ind_nbr(nn)

             ! list of unique grid cells, starting with the grid cell that contains lat/lon

             N_cell       = 1

             i_cell_c3(1) = i_ind(ind_nbr(nn))
             j_cell_c3(1) = j_ind(ind_nbr(nn))

             do k=1,N_pts_nbr(nn)

                kk = kk + 1

                if (any( i_cell_c3(1:N_cell)==i_ind_nbr(kk) .and.               &
                         j_cell_c3(1:N_cell)==j_ind_nbr(kk)       ))  cycle  ! already in list

                N_cell            = N_cell + 1

                i_cell_c3(N_cell) = i_ind_nbr(kk)
                j_cell_c3(N_cell) = j_ind_nbr(kk)

             end do

             ! Step 2.): tiles within additional grid cells (grid cell 1 was searched in Step 1.)

             min_dist = 1.e10    ! initialize

             if (N_cell>1)                                                       &
                  call get_tile_num_in_bbox_c3( lat(n), lon(n), N_cell-1,        &
                  i_cell_c3(2:N_cell), j_cell_c3(2:N_cell),                      &
                  tile_coord, tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij, &
                  tile_num(n), min_dist, min_dist_x(ind_nbr(nn)), min_dist_y(ind_nbr(nn)) )

             ! Step 3.): no bounding box contains lat/lon:  find closest tile among those
             !           in all grid cells identified above that are within max_dist

             if (tile_num(n)==nodata_tilenum) then

                do mm=1,N_cell

                   do k=1,N_tile_in_cell_ij(i_cell_c3(mm),j_cell_c3(mm))

                      this_tile_num = tile_num_in_cell_ij(i_cell_c3(mm),j_cell_c3(mm),k)

                      tmp_dist_x = abs(lon(n) - tile_coord(this_tile_num)%com_lon)
                      tmp_dist_y = abs(lat(n) - tile_coord(this_tile_num)%com_lat)

                      ! skip tile that is too far away (see bounding box and max_dist check below)

                      if (tmp_dist_x > max_dist_x(n) .or. tmp_dist_y > max_dist_y)  cycle

                      tmp_dist   = tmp_dist_x + tmp_dist_y    ! "closest": Minkowski norm

                      if (tmp_dist<min_dist) then

                         min_dist                 = tmp_dist
                         min_dist_x(ind_nbr(nn))  = tmp_dist_x
                         min_dist_y(ind_nbr(nn))  = tmp_dist_y

                         tile_num(n)              = this_tile_num

                      end if

                   end do

                end do

             end if

          end do

          nn_start = nn_end + 1

       end do

       ! check whether lat/lon is inside bounding box of given tile
       ! (see get_tile_num_from_latlon())

       do nn=1,N_tmp

          n = ind_start - 1 + nn

          if (tile_num(n)>0) then

             outside_bbox  = (                                          &
                  lon(n) < tile_coord(tile_num(n))%min_lon    .or.      &
                  lon(n) > tile_coord(tile_num(n))%max_lon    .or.      &
                  lat(n) < tile_coord(tile_num(n))%min_lat    .or.      &
                  lat(n) > tile_coord(tile_num(n))%max_lat           )

             too_far_away = (                                           &
                  min_dist_x(nn) > max_dist_x(n)              .or.      &
                  min_dist_y(nn) > max_dist_y                        )

             ! keep tile_num unless obs is outside the bounding box *and* too far away

             if (outside_bbox .and. too_far_away)  tile_num(n) = nodata_tilenum

          end if

       end do

    end do

    deallocate( i_ind, j_ind, ind_nbr, N_pts_nbr, min_dist_x, min_dist_y )
    deallocate( i_ind_nbr, j_ind_nbr, tmp_lat, tmp_lon )

  end subroutine get_tile_num_from_latlon_c3

  ! *******************************************************************

  subroutine get_tile_num_in_bbox_c3( lat, lon, N_cell, i_cell, j_cell,        &
       tile_coord, tile_grid, N_tile_in_cell_ij, tile_num_in_cell_ij,          &
       tile_num, min_dist, min_dist_x, min_dist_y )

    ! among the tiles within the given grid cells of the cube-sphere tile_grid, find the
    !  tile whose bounding box contains lat/lon and that has the minimum normalized square
    !  distance (dx/hx)**2 + (dy/hy)**2 (< min_dist), where hx and hy are half the width
    !  and height of the bounding box; see get_tile_num_from_latlon_c3()
    !
    ! tile_num, min_dist, min_dist_x, and min_dist_y change only if such a tile is found
    !
    ! amfox, 8 Nov 2026 (moved here from get_tile_num_from_latlon())

    implicit none

    real,                                                intent(in)    :: lat, lon

    integer,                                             intent(in)    :: N_cell

    integer, dimension(N_cell),                          intent(in)    :: i_cell, j_cell

    type(tile_coord_type), dimension(:),                 pointer       :: tile_coord ! input

    type(grid_def_type),                                 intent(in)    :: tile_grid

    integer, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in)    :: N_tile_in_cell_ij

    integer, dimension(:,:,:),                           pointer       :: tile_num_in_cell_ij ! input

    integer,                                             intent(inout) :: tile_num

    real,                                                intent(inout) :: min_dist, min_dist_x, min_dist_y

    ! local variables

    integer                      :: k, kk, this_tile_num

    real                         :: tmp_dist, tmp_dist_x, tmp_dist_y, tmp_hx, tmp_hy

    ! -----------------------------------------------------------

    do kk=1,N_cell

       do k=1,N_tile_in_cell_ij(i_cell(kk),j_cell(kk))

          this_tile_num = tile_num_in_cell_ij(i_cell(kk),j_cell(kk),k)

          if ( tile_coord(this_tile_num)%min_lon <= lon    .and.      &
               tile_coord(this_tile_num)%max_lon >= lon    .and.      &
               tile_coord(this_tile_num)%min_lat <= lat    .and.      &
               tile_coord(this_tile_num)%max_lat >= lat           ) then

             ! lat/lon is inside bounding box of tile

             tmp_hx     = 0.5*(tile_coord(this_tile_num)%max_lon - tile_coord(this_tile_num)%min_lon)
             tmp_hy     = 0.5*(tile_coord(this_tile_num)%max_lat - tile_coord(this_tile_num)%min_lat)

             tmp_hx     = max( tmp_hx, 1.e-6 )   ! avoid division by zero
             tmp_hy     = max( tmp_hy, 1.e-6 )

             tmp_dist_x = abs(lon - tile_coord(this_tile_num)%com_lon)
             tmp_dist_y = abs(lat - tile_coord(this_tile_num)%com_lat)

             tmp_dist   = (tmp_dist_x/tmp_hx)**2 + (tmp_dist_y/tmp_hy)**2

             if (tmp_dist<min_dist) then

                min_dist    = tmp_dist
                min_dist_x  = tmp_dist_x
                min_dist_y  = tmp_dist_y

                tile_num    = this_tile_num

             end if

          end if

       end do

    end do

  end subroutine get_tile_num_in_bbox_c3

  ! *******************************************************************

  subroutine get_ij_ind_from_latlon( tile_grid, lat, lon, i_ind, j_ind )

    ! NOTE order of input arguments ("lat", "lon", "lon_ind", "lat_ind")

    ! find i/j_ind of grid cell that contains lat/lon
    !
    ! that is, for given lat/lon, find out corresponding i/j_ind into an array
    ! of size N_lon-by-N_lat whose lat/lon coordinates are defined by tile_grid
    
    ! NOTE: ALL arrays in LDASsa are declared with base 1
    
    ! typical use: have lat/lon of an obs, want to know i/j_ind within the 
    ! grid (array) that defines the tile space ("tile_grid_d")

    ! major revision by reichle, 11 May 2011
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    real, intent(in) :: lat, lon
    integer, intent(out) :: i_ind, j_ind
    
    real    :: lats(1),    lons(1)
    integer :: i_inds(1), j_inds(1)

    ! local variables
    
    real :: r, s, i_indg, j_indg

    character(len=*), parameter :: Iam = 'get_ij_ind_from_latlon'
    character(len=400) :: err_msg

    ! -------------------------------------------------------------
    !
    ! select grid type
    
    if     (index(tile_grid%gridtype, 'EASEv')/=0) then

       ! EASE grid lat/lon to index provides *global*, *0-based* index!

       call MAPL_ease_convert(tile_grid%gridtype, lat, lon, r, s)
       
       i_indg = nint(r)    ! i_ind or lon_ind
       j_indg = nint(s)    ! j_ind or lat_ind

       ! convert to index into array defined by tile_grid_d
       
       i_ind = i_indg - tile_grid%i_offg - (tile_grid%ind_base - 1)
       j_ind = j_indg - tile_grid%j_offg - (tile_grid%ind_base - 1)
       
    elseif (index(tile_grid%gridtype, 'LatLon')/=0) then
       
       ! ll_lon and ll_lat refer to lower left  corner of grid cell
       ! ur_lon and ur_lat refer to upper right corner of grid cell
       ! (as opposed to the grid point in the center of the grid cell)
       !
       ! ALL arrays in LDASsa are declared with base 1 (--> use "ceiling")

       if (tile_grid%i_dir==1) then
          
          i_ind = ceiling( (lon - tile_grid%ll_lon)/tile_grid%dlon )
          
       else
          
          i_ind = ceiling( (tile_grid%ur_lon - lon)/tile_grid%dlon )

       end if


       if (tile_grid%j_dir==1) then
          
          j_ind = ceiling( (lat - tile_grid%ll_lat)/tile_grid%dlat )
          
       else
          
          j_ind = ceiling( (tile_grid%ur_lat - lat)/tile_grid%dlat )
          
       end if

    elseif (index(tile_grid%gridtype, 'c3')/=0) then

       ! cube-sphere grid requires the ESMF grid (amfox, 8 Nov 2026)

       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'cube-sphere grid: use get_ij_ind_from_latlon_c3()')

    else

       call ldas_abort(LDAS_GENERIC_ERROR, Iam, 'unknown grid type')

    end if
    
  end subroutine get_ij_ind_from_latlon
  
  ! *******************************************************************

  subroutine get_ij_ind_from_latlon_c3( tile_grid_g_ESMF, tile_grid, N_latlon, lat, lon, i_ind, j_ind )

    ! find i/j_ind of cube-sphere ("c3") grid cells that contain lat/lon
    !
    ! as get_ij_ind_from_latlon(), but for a vector of lat/lon and cube-sphere tile_grid;
    !  returns i/j_ind into an array of size N_lon-by-N_lat defined by tile_grid,
    !  limited to 1,...,N_lon and 1,...,N_lat
    !
    ! MAPL_GetGlobalHorzIJIndex() needs the ESMF grid and returns *global*, *1-based*
    !  indices i=1,...,N and j=1,...,6*N (6 faces stacked in j), as in tile_coord%i_indg
    !  and tile_coord%j_indg.  It fails for grids that it cannot handle.  (It also handles
    !  stretched grids, but init_tile_grid_g_ESMF_c3() rejects them, see there.)
    !
    ! NOTE: Each call of MAPL_GetGlobalHorzIJIndex() checks the grid, and the MAPL subroutine
    !       allocates ~100 bytes per lat/lon --> call with many, but not too many, lat/lon.
    !
    ! NOTE: Before this subroutine is called by root_proc only (e.g., in obs readers),
    !       init_tile_grid_g_ESMF_c3() must have been called by all processors.
    !
    ! amfox, 8 Nov 2026

    implicit none

    type(ESMF_Grid),              intent(in)  :: tile_grid_g_ESMF

    type(grid_def_type),          intent(in)  :: tile_grid

    integer,                      intent(in)  :: N_latlon

    real,    dimension(N_latlon), intent(in)  :: lat, lon      ! [deg]

    integer, dimension(N_latlon), intent(out) :: i_ind, j_ind

    ! local variables

    type(ESMF_Grid)                       :: tmp_grid

    real*8,  dimension(:),   allocatable  :: tmp_lat, tmp_lon

    integer                               :: status

    character(len=*),  parameter :: Iam = 'get_ij_ind_from_latlon_c3'
    character(len=400)           :: err_msg

    ! -------------------------------------------------------------

    if (N_latlon<1) return

    ! convert lat/lon from [deg] to [rad]

    allocate( tmp_lat(N_latlon), tmp_lon(N_latlon) )

    tmp_lat = real(lat,8)*(MAPL_PI_R8/180.d0)
    tmp_lon = real(lon,8)*(MAPL_PI_R8/180.d0)

    ! MAPL_GetGlobalHorzIJIndex() needs intent(inout) for the grid (to store the grid corners
    !  as attributes, see init_tile_grid_g_ESMF_c3()) --> pass a copy of the ESMF grid handle (same grid)

    tmp_grid = tile_grid_g_ESMF

    call MAPL_GetGlobalHorzIJIndex( N_latlon, i_ind, j_ind, lonR8=tmp_lon, latR8=tmp_lat,   &
         Grid=tmp_grid, rc=status )

    if (status/=0) then

       err_msg = 'MAPL_GetGlobalHorzIJIndex() failed'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if

    ! convert to index into array defined by tile_grid

    i_ind = i_ind - tile_grid%i_offg - (tile_grid%ind_base - 1)
    j_ind = j_ind - tile_grid%j_offg - (tile_grid%ind_base - 1)

    ! make sure that i/j_ind is within bounds

    i_ind = min( max(i_ind, 1), tile_grid%N_lon )
    j_ind = min( max(j_ind, 1), tile_grid%N_lat )

    deallocate( tmp_lat, tmp_lon )

  end subroutine get_ij_ind_from_latlon_c3

  ! *******************************************************************

  subroutine init_tile_grid_g_ESMF_c3( tile_grid_g_ESMF, tile_grid )

    ! prepare cube-sphere ("c3") ESMF grid for get_ij_ind_from_latlon_c3() and check that
    !  the grid matches the global cube-sphere tile_grid (tile_grid_g)
    !
    ! MUST be called by ALL processors (once), before get_ij_ind_from_latlon_c3() is
    !  called by root_proc only (e.g., in obs readers):
    !
    ! MAPL_GetGlobalHorzIJIndex() checks the grid with MAPL_GridGetCorners().  The first
    !  call of MAPL_GridGetCorners() for a given grid needs all processors (halo exchange)
    !  and stores the grid corners as attributes of the grid.  Later calls use these
    !  attributes and are local.  Here, MAPL_GetGlobalHorzIJIndex() is called for zero
    !  lat/lon, which only checks the grid.
    !
    ! Stretched cube-sphere grids are not supported:  MAPL_GetGlobalHorzIJIndex() handles
    !  them, but the search for nearby grid cells in get_tile_num_from_latlon_c3() uses the
    !  grid spacing of an unstretched grid (tile_grid%dlon) and misses many tiles where the
    !  stretched grid is finer.  MAPL marks a stretched grid with the attributes
    !  STRETCH_FACTOR, TARGET_LON and TARGET_LAT (see MAPL_Reverse_Schmidt()).
    !
    ! amfox, 8 Nov 2026

    implicit none

    type(ESMF_Grid),              intent(in)  :: tile_grid_g_ESMF

    type(grid_def_type),          intent(in)  :: tile_grid

    ! local variables

    type(ESMF_Grid)                       :: tmp_grid

    real*8,  dimension(0)                 :: tmp_lat, tmp_lon

    integer, dimension(0)                 :: tmp_i_ind, tmp_j_ind

    integer, dimension(3)                 :: dims

    integer                               :: status

    logical                               :: has_factor, has_lon, has_lat

    character(len=*),  parameter :: Iam = 'init_tile_grid_g_ESMF_c3'
    character(len=400)           :: err_msg

    ! -------------------------------------------------------------

    tmp_grid = tile_grid_g_ESMF         ! copy of ESMF grid handle (same grid), see get_ij_ind_from_latlon_c3()

    ! check grid size (N-by-6*N for cube-sphere grid)

    dims = 0

    call MAPL_GridGet( tmp_grid, globalCellCountPerDim=dims, rc=status )

    if (status/=0) then

       err_msg = 'MAPL_GridGet() failed'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if

    if ( (dims(1)/=tile_grid%N_lon) .or. (dims(2)/=tile_grid%N_lat) ) then

       write (err_msg,'(A,2i8,A,2i8)') 'size of ESMF grid ', dims(1:2),       &
            ' does not match size of cube-sphere tile_grid ', tile_grid%N_lon, tile_grid%N_lat
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if

    ! reject stretched grid (see above)

    call ESMF_AttributeGet( tmp_grid, name='STRETCH_FACTOR', isPresent=has_factor, rc=status )
    if (status==0)                                                                       &
         call ESMF_AttributeGet( tmp_grid, name='TARGET_LON', isPresent=has_lon,    rc=status )
    if (status==0)                                                                       &
         call ESMF_AttributeGet( tmp_grid, name='TARGET_LAT', isPresent=has_lat,    rc=status )

    if (status/=0) then

       err_msg = 'ESMF_AttributeGet() failed'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if

    if (has_factor .and. has_lon .and. has_lat) then

       err_msg = 'stretched cube-sphere grid not supported for mapping obs to tiles (obs assimilation)'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if

    ! check grid and store grid corners as grid attributes (on all processors)

    call MAPL_GetGlobalHorzIJIndex( 0, tmp_i_ind, tmp_j_ind, lonR8=tmp_lon, latR8=tmp_lat,  &
         Grid=tmp_grid, rc=status )

    if (status/=0) then

       err_msg = 'MAPL_GetGlobalHorzIJIndex() failed; cube-sphere grid not supported'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)

    end if

    if (logit) write (logunit,*) trim(Iam), ': cube-sphere tile_grid_g_ESMF is ready for obs-to-tile mapping, size = ', dims(1:2)
    if (logit) write (logunit,*)

  end subroutine init_tile_grid_g_ESMF_c3

  ! *******************************************************************

  subroutine get_tile_num_in_ellipse( lon, lat, r_x, r_y,                 &
       N_tile, tile_coord, tile_grid,                                     &
       N_tile_in_cell_ij, tile_num_in_cell_ij,                            &
       N_tile_in_ellipse,  tile_num_in_ellipse,  norm_square_distance   )
            
    ! reichle, 22 Feb 2015

    ! return tile numbers of all tiles within ellipse with center at lat/lon
    !  and axes r_x, r_y
    !
    ! also return normalized square distance of tile center-of-mass lat/lon from lat/lon

    ! TO DO: what about dateline? poles?  [ignore for now]
    
    ! NOTE: get_tile_num_in_ellipse() returns zero tiles (N_tile_in_ellipse=0)
    !        if ellipse straddles date line and EASE[v2] grid is used
    !       - reichle, 17 Apr 2017
    
    
    implicit none
    
    real,                                     intent(in)    :: lon, lat, r_x, r_y

    integer,                                  intent(in)    :: N_tile
    
    type(tile_coord_type), dimension(:),      pointer       :: tile_coord ! input
    
    type(grid_def_type),                      intent(in)    :: tile_grid
    
    integer,               dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in) :: &
         N_tile_in_cell_ij
        
    integer,               dimension(:,:,:),  pointer       :: tile_num_in_cell_ij ! input
    
    integer,                                  intent(out)   :: N_tile_in_ellipse
        
    integer,               dimension(N_tile), intent(out)   :: tile_num_in_ellipse

    real,                  dimension(N_tile), intent(out)   :: norm_square_distance
    
    ! local variables
    
    real    :: minlat, maxlat, minlon, maxlon, tmp_dx, tmp_dy, tmp_dist
    real    :: r_x_square, r_y_square
    
    integer :: i_ind_ll, i_ind_ur, i_min, i_max
    integer :: j_ind_ll, j_ind_ur, j_min, j_max

    integer :: ii, jj, kk, mm, nn
    
    ! -------------------------------------------------------------
    !
    ! initialize
    
    N_tile_in_ellipse   = 0
    
    tile_num_in_ellipse = nodata_generic
    
    r_x_square          = r_x**2
    r_y_square          = r_y**2
    
    ! identify rectangle +/- r_x and +/- r_y around lat/lon
    
    tmp_dx = tile_grid%dlon  ! broaden rectangle to avoid missing tiles
    tmp_dy = tile_grid%dlon  ! use *dlon* because for EASE grid *dlat* is ill-defined   
    
    minlon = lon - r_x - tmp_dx
    maxlon = lon + r_x + tmp_dx
    minlat = lat - r_y - tmp_dy
    maxlat = lat + r_y + tmp_dy
    
    ! find i,j indices of two opposite corners (lower left and upper right)
    
    call get_ij_ind_from_latlon( tile_grid, minlat, minlon, i_ind_ll, j_ind_ll )
    call get_ij_ind_from_latlon( tile_grid, maxlat, maxlon, i_ind_ur, j_ind_ur )
    
    ! restrict i_ind_* and j_ind_* to fall within tile_grid
    
    i_ind_ll = min( tile_grid%N_lon, max( 1, i_ind_ll ))
    j_ind_ll = min( tile_grid%N_lat, max( 1, j_ind_ll ))
    i_ind_ur = min( tile_grid%N_lon, max( 1, i_ind_ur ))
    j_ind_ur = min( tile_grid%N_lat, max( 1, j_ind_ur ))

    ! find smaller and larger of the two indices
    ! (index may run north-to-south/west-to-east depending on i_dir/j_dir)
    
    i_min = min( i_ind_ll, i_ind_ur )
    i_max = max( i_ind_ll, i_ind_ur )

    j_min = min( j_ind_ll, j_ind_ur )
    j_max = max( j_ind_ll, j_ind_ur )
    
    ! loop through all grid cells "inside" rectangle and identify tiles w/in ellipse

    norm_square_distance = nodata_generic

    mm = 0
    
    do ii=i_min,i_max
       do jj=j_min,j_max
          
          do kk=1,N_tile_in_cell_ij(ii,jj)
             
             nn = tile_num_in_cell_ij(ii,jj,kk)
             
             tmp_dist = (                                            &
                  (tile_coord(nn)%com_lon - lon)**2 / r_x_square +    &
                  (tile_coord(nn)%com_lat - lat)**2 / r_y_square   )
             
             if (tmp_dist <= 1.) then
                
                ! record tile number and distance
                
                mm = mm+1
                
                tile_num_in_ellipse(mm) = nn
                
                norm_square_distance(mm) = tmp_dist
                
             end if
                          
          end do
          
       end do
    end do
    
    N_tile_in_ellipse = mm
    
  end subroutine get_tile_num_in_ellipse

  ! **********************************************************************
  
  subroutine tile2grid_simple( N_tile, i_indgs, j_indgs, tile_grid, tile_data, grid_data)

    ! Interpolate from tile space to grid space without interpolation/weighted/no-data-handling.
    ! Simply assign the tile value to the grid cell (last assignment prevails)

    implicit none
    
    integer, intent(in) :: N_tile

    integer, intent(in) :: i_indgs(:)     ! dimension(N_tile)
    integer, intent(in) :: j_indgs(:)     ! dimension(N_tile)
     
    type(grid_def_type), intent(in) :: tile_grid
    
    real, dimension(N_tile), intent(in) :: tile_data
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(out) :: grid_data
        
    ! local variables
    
    integer :: n, i, j, off_i, off_j
    real, parameter :: no_data= -9999.   
    character(len=*), parameter :: Iam = 'tile2grid_simple'
    character(len=400) :: err_msg

    ! ------------------------------------
    
    if (size(i_indgs)/=N_tile .or. size(j_indgs)/=N_tile) then
       err_msg = '[i,j]_indg and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    !
    ! adjust for 0-based indexing (eg., EASE grids)
    !
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    ! loop through tile space

    grid_data = no_data
    do n=1,N_tile
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j
       grid_data(i,j) = tile_data(n)
    end do
    
  end subroutine tile2grid_simple

  ! *******************************************************************
  !
  ! The subroutine tile2grid_full() is the original LDASsa tile2grid() 
  !  subroutine with weighted averaging and no-data-handling.
  !
  ! This subroutine should no longer be used and is provided here for
  !   reference only.  Use MAPL LocationStream instead.
  !
  ! - reichle, 13 July 2020

  subroutine tile2grid_full( N_tile, tile_coord, tile_grid, tile_data,    &
       grid_data, no_data_value, no_data_tol, echo )
    
    ! map from tile space to tile definition grid
    !
    ! NOTE: tile_coord must match tile_data
    !
    ! optional inputs: 
    !   no_data_value :
    !   no_data_tol   : tolerance when checking tile_data 
    !                     against no_data_value
    !   echo          : echo no_data_value and tolerance
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domaim of interest.  
    ! With these offsets tile2grid() can be used to map from tile space to a 
    ! subgrid of "tile_grid_g"
    
    ! reichle, 28 Jan 2005
    ! reichle, 16 May 2005 - added offset for "domain" grid
    ! reichle, 21 May 2010 - off_i, off_j now part of grid_def_type

    implicit none
    
    integer, intent(in) :: N_tile
    
    type(tile_coord_type), dimension(:), pointer :: tile_coord  ! input
    
    type(grid_def_type), intent(in) :: tile_grid
    
    real, dimension(N_tile), intent(in) :: tile_data
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(out) :: grid_data
    
    real, intent(in), optional :: no_data_value, no_data_tol
    
    logical, intent(in), optional :: echo
        
    ! local variables
    
    integer :: n, i, j, off_i, off_j
    
    real :: w, no_data, tol
    
    real, parameter :: no_data_default     = -9999.
    real, parameter :: no_data_tol_default = 1e-4
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat) :: wgrid
    
    character(len=*), parameter :: Iam = 'tile2grid'
    character(len=400) :: err_msg

    ! ------------------------------------
    
    if (size(tile_coord)/=N_tile) then
       err_msg = 'tile_coord and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    if (present(no_data_value)) then
       no_data = no_data_value
    else
       no_data = no_data_default
    end if
    
    if (present(no_data_tol)) then
       tol = no_data_tol
    else
       tol = no_data_tol_default
    end if
    
    if (present(echo)) then
       if (echo .and. logit) then
          write (logunit,*) 'tile2grid: using no-data-value = ', no_data , &
               ' with tolerance = ', tol
       end if
    end if
    
    ! ------------------------------------------------------
    !
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
    
    ! initialize
    
    grid_data = 0.
    wgrid     = 0.
    
    ! loop through tile space
    
    do n=1,N_tile
       
       i = tile_coord(n)%i_indg - off_i
       j = tile_coord(n)%j_indg - off_j
       
       w = tile_coord(n)%frac_cell
       
       if (abs(tile_data(n)-no_data)>tol) then
          
          grid_data(i,j) = grid_data(i,j) + w*tile_data(n)
          
          wgrid(i,j) = wgrid(i,j) + w
          
       end if
       
    end do
    
    ! normalize and set no-data-value
    
    do i=1,tile_grid%N_lon
       do j=1,tile_grid%N_lat
          
          if (wgrid(i,j)>0.) then
             
             grid_data(i,j) = grid_data(i,j)/wgrid(i,j)
             
          else
             
             grid_data(i,j) = no_data
             
          end if
          
       end do
    end do
    
  end subroutine tile2grid_full
  
  ! *******************************************************************

  subroutine tile_mask_grid( tile_grid, N_tile, i_indgs,j_indgs, grid_data)
    
    ! set grid cell to no value if there is no tile in it
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    integer, intent(in) :: N_tile
    
    integer, intent(in) :: i_indgs(:)
    integer, intent(in) :: j_indgs(:)
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(inout) :: grid_data
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat) :: grid
        
    ! local variables
    real, parameter :: no_data= -9999.
    integer :: n, i, j, off_i, off_j

    character(len=*), parameter :: Iam = 'tile_mask_grid'
    character(len=400) :: err_msg
    
    ! ------------------------------------
    
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
       
    grid = no_data
    do n=1,N_tile
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j
       grid(i,j) = grid_data(i,j)
    end do
    grid_data = grid

  end subroutine tile_mask_grid

  ! **********************************************************************

  subroutine grid2tile_real( tile_grid, N_tile, i_indgs,j_indgs, grid_data, tile_data)
    
    ! map from grid to tile space
    !
    ! The indices tile_coord%i_indg and tile_coord%j_indg refer to the *global*
    ! tile definition grid (as obtained from the tile_coord_file).
    ! Integers "off_i" and "off_j" describe the offset between the global 
    ! "tile_grid_g" and a smaller "tile_grid_d" for the domain of interest.  
    ! With these offsets grid2tile() can be used to map from a
    ! subgrid of "tile_grid_g" to tile space
    !
    ! reichle, 28 Jan 2005
    ! reichle, 16 May 2005 - added offset for "domain" grid
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    integer, intent(in) :: N_tile
    
    integer, intent(in) :: i_indgs(:)
    integer, intent(in) :: j_indgs(:)
    
    real, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in) :: grid_data
    
    real, dimension(N_tile), intent(out) :: tile_data
        
    ! local variables
    
    integer :: n, i, j, off_i, off_j

    character(len=*), parameter :: Iam = 'grid2tile_real'
    character(len=400) :: err_msg
    
    ! ------------------------------------
    
    if (size(i_indgs)/=N_tile .or. size(j_indgs)/=N_tile) then
       err_msg = '[i,j]_indg and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
    
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
        
    do n=1,N_tile
       
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j

       tile_data(n) = grid_data(i,j)
       
    end do
    
  end subroutine grid2tile_real

  ! **********************************************************************
  
  subroutine grid2tile_real8( tile_grid, N_tile, i_indgs,j_indgs, grid_data_8, tile_data_8)
    
    ! same as grid2tile_real but for real*8
    !
    ! reichle,  3 Feb 2014
    
    implicit none
    
    type(grid_def_type), intent(in) :: tile_grid
    
    integer, intent(in) :: N_tile
    integer, intent(in) :: i_indgs(:)
    integer, intent(in) :: j_indgs(:)
    
    real*8, dimension(tile_grid%N_lon,tile_grid%N_lat), intent(in) :: grid_data_8
    
    real*8, dimension(N_tile), intent(out) :: tile_data_8
    
    ! local variables
    
    integer :: n, i, j, off_i, off_j
    
    character(len=*), parameter :: Iam = 'grid2tile_real8'
    character(len=400) :: err_msg

    ! ------------------------------------
    
    if (size(i_indgs)/=N_tile .or. size(j_indgs)/=N_tile) then
       err_msg = '[i,j]_indg and tile_data do not match.'
       call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if
        
    ! adjust for 0-based indexing (eg., EASE grids)
    
    off_i = tile_grid%i_offg + (tile_grid%ind_base - 1)
    off_j = tile_grid%j_offg + (tile_grid%ind_base - 1)
        
    do n=1,N_tile
       
       i = i_indgs(n) - off_i
       j = j_indgs(n) - off_j
       
       tile_data_8(n) = grid_data_8(i,j)
       
    end do
    
  end subroutine grid2tile_real8

  ! **********************************************************************

end module LDAS_TileCoordRoutines


! ================================= EOF =======================================
