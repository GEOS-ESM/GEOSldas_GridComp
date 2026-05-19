#include "MAPL_Generic.h"

module GEOS_RouteAvgGridCompMod
  use ESMF
  use MAPL_Mod
  use GEOS_RouteGridCompMod, only : pfaf_grid, pfaf_locstream

  implicit none

  private

  public   :: SetServices

  integer  :: NUM_ENSEMBLE
  integer  :: collect_route_counter
  
contains
  
  subroutine SetServices(gc, rc)
    type(ESMF_GridComp), intent(inout) :: gc ! gridded component
    integer, optional                  :: rc ! return code

    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! Get my name and setup traceback handle
    Iam = 'SetServices'
    call ESMF_GridCompGet(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::" // Iam

    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_INITIALIZE,                                                &
         Initialize,                                                            &
         rc=status                                                              &
         )
    VERIFY_(status)

    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_RUN,                                                       &
         RUN,                                                                   &
         rc=status                                                              &
         )
    VERIFY_(status)

    call MAPL_GridCompSetEntryPoint(                                            &
         gc,                                                                    &
         ESMF_METHOD_FINALIZE,                                                  &
         Finalize,                                                              &
         rc=status                                                              &
         )
    VERIFY_(status)

!  export route

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'total_runoff_flow_from_catchment'                 ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QRUNOFF_CAT'              ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'discharge_from_local_streams_to_main_river'       ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QSFLOW'                   ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'discharge_from_upstream_catchments_to_main_river' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QINFLOW'                  ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'discharge_from_main_river'                        ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QOUTFLOW'                 ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'discharge_from_reservoir'                         ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QRES'                     ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'discharge_to_downstream_catchment'                ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QCAT'                     ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )

    ! Set profiling timers
    call MAPL_TimerAdd(gc, name="Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Collect_route", rc=status)
    VERIFY_(status)

    ! Call SetServices for children
    call MAPL_GenericSetServices(gc, rc=status)
    VERIFY_(status)
    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

  subroutine Initialize(gc, import, export, clock, rc)
    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null()

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Initialize"

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

   ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Initialize")

    call MAPL_GetResource ( MAPL, NUM_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)

    collect_route_counter    = 0

    call MAPL%grid%set(pfaf_grid, _RC)
    call ESMF_GridCompSet(gc, grid=pfaf_grid, RC=status)
    VERIFY_(STATUS)
    call MAPL_set(MAPL, locstream = pfaf_locstream, rc=status)
    VERIFY_(STATUS)

    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    call MAPL_TimerOff(MAPL, "Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")  

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

  subroutine RUN (gc, import, export, clock, rc)
    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! MAPL variables
    type(MAPL_MetaComp), pointer :: MAPL=>null() ! MAPL obj

    ! Pointers to imports
    real, pointer :: QSFLOW(:)=>null()
    real, pointer :: QOUTFLOW(:)=>null()
    real, pointer :: QRES(:)=>null()
    real, pointer :: QCAT(:)=>null()
    real, pointer :: QINFLOW(:)=>null()
    real, pointer :: QRUNOFF_CAT(:)=>null()

    ! Pointers to exports
    real, pointer :: QSFLOW_enavg(:)=>null()
    real, pointer :: QOUTFLOW_enavg(:)=>null()
    real, pointer :: QRES_enavg(:)=>null()
    real, pointer :: QCAT_enavg(:)=>null()
    real, pointer :: QINFLOW_enavg(:)=>null()
    real, pointer :: QRUNOFF_CAT_enavg(:)=>null()

     ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_route_ens_averaging"

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Collect_route")

    call MAPL_GetPointer(import, QSFLOW,          'QSFLOW' ,     _RC)
    call MAPL_GetPointer(import, QOUTFLOW,        'QOUTFLOW',    _RC)
    call MAPL_GetPointer(import, QRES,            'QRES',        _RC)
    call MAPL_GetPointer(import, QCAT,            'QCAT',        _RC)
    call MAPL_GetPointer(import, QINFLOW,         'QINFLOW',     _RC)
    call MAPL_GetPointer(import, QRUNOFF_CAT,     'QRUNOFF_CAT', _RC)

    call MAPL_GetPointer(export, QSFLOW_enavg,          'QSFLOW' ,     _RC)
    call MAPL_GetPointer(export, QOUTFLOW_enavg,        'QOUTFLOW',    _RC)
    call MAPL_GetPointer(export, QRES_enavg,            'QRES',        _RC)
    call MAPL_GetPointer(export, QCAT_enavg,            'QCAT',        _RC)
    call MAPL_GetPointer(export, QINFLOW_enavg,         'QINFLOW',     _RC)
    call MAPL_GetPointer(export, QRUNOFF_CAT_enavg,     'QRUNOFF_CAT', _RC)

    if(collect_route_counter ==0) then
       if(associated(QSFLOW_enavg))       QSFLOW_enavg      = 0.0
       if(associated(QOUTFLOW_enavg))     QOUTFLOW_enavg    = 0.0
       if(associated(QRES_enavg))         QRES_enavg        = 0.0
       if(associated(QCAT_enavg))         QCAT_enavg        = 0.0
       if(associated(QINFLOW_enavg))      QINFLOW_enavg     = 0.0
       if(associated(QRUNOFF_CAT_enavg))  QRUNOFF_CAT_enavg = 0.0
    endif

    if(associated(QSFLOW_enavg))       QSFLOW_enavg      = QSFLOW_enavg      + QSFLOW
    if(associated(QOUTFLOW_enavg))     QOUTFLOW_enavg    = QOUTFLOW_enavg    + QOUTFLOW
    if(associated(QRES_enavg))         QRES_enavg        = QRES_enavg        + QRES
    if(associated(QCAT_enavg))         QCAT_enavg        = QCAT_enavg        + QCAT
    if(associated(QINFLOW_enavg))      QINFLOW_enavg     = QINFLOW_enavg     + QINFLOW
    if(associated(QRUNOFF_CAT_enavg))  QRUNOFF_CAT_enavg = QRUNOFF_CAT_enavg + QRUNOFF_CAT

    collect_route_counter = collect_route_counter + 1

    if (collect_route_counter == NUM_ENSEMBLE ) then
       collect_route_counter = 0
       if(associated(QSFLOW_enavg))       QSFLOW_enavg      = QSFLOW_enavg  /NUM_ENSEMBLE
       if(associated(QOUTFLOW_enavg))     QOUTFLOW_enavg    = QOUTFLOW_enavg/NUM_ENSEMBLE
       if(associated(QRES_enavg))         QRES_enavg        = QRES_enavg    /NUM_ENSEMBLE
       if(associated(QCAT_enavg))         QCAT_enavg        = QCAT_enavg    /NUM_ENSEMBLE
       if(associated(QINFLOW_enavg))      QINFLOW_enavg     = QINFLOW_enavg /NUM_ENSEMBLE
       if(associated(QRUNOFF_CAT_enavg))  QRUNOFF_CAT_enavg = QRUNOFF_CAT_enavg/NUM_ENSEMBLE
    endif

    call MAPL_TimerOff(MAPL, "Collect_route")
    call MAPL_TimerOff(MAPL, "TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

  subroutine Finalize(gc, import, export, clock, rc)
    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name

    ! Local variables
    type(MAPL_MetaComp), pointer :: MAPL=>null() ! MAPL obj

    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Finalize"

   ! Call Finalize for every child
    call MAPL_GenericFinalize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Finalize

end module GEOS_RouteAvgGridCompMod
