#include "MAPL_Generic.h"

module GEOS_LakeAvgGridCompMod

  use ESMF
  use MAPL_Mod
  
  implicit none
  
  private

  public   :: SetServices
  
  integer  :: NUM_ENSEMBLE_LAKE = 1
  integer  :: collect_lake_counter
  
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

    ! subset of exports from Lake GridComp, with *updated* long_name attributes
    
    call MAPL_AddExportSpec(GC,                            &
         LONG_NAME          = 'total_evaporation_lake'    ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'EVAPOUT'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                            &
         LONG_NAME          = 'runoff_total_flux_lake'    ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RUNOFF'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                            &
         LONG_NAME          = 'sensible_heat_flux_lake'   ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'SHOUT'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC                            ,&
         LONG_NAME          = 'net_longwave_flux_lake'    ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'LWNDSRF'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC                            ,&
         LONG_NAME          = 'net_shortwave_flux_lake'   ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'SWNDSRF'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                            &
         LONG_NAME          = 'latent_heat_flux_lake'     ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'HLATN'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         LONG_NAME          = 'surface_temperature_lake',          &
         UNITS              = 'K',                                 &
         SHORT_NAME         = 'TST',                               &
         DIMS               = MAPL_DimsTileOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         LONG_NAME          = 'surface_specific_humidity_lake',    &
         UNITS              = 'kg kg-1',                           &
         SHORT_NAME         = 'QST',                               &
         DIMS               = MAPL_DimsTileOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

  
    ! Set profiling timers
    call MAPL_TimerAdd(gc, name="Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Collect_lake", rc=status)
    VERIFY_(status)

    ! Call SetServices for children
    call MAPL_GenericSetServices(gc, rc=status)
    VERIFY_(status)
    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

  ! -------------------------------------------------------------------------------------------
  
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

    !call MAPL_GetResource ( MAPL, NUM_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    !VERIFY_(STATUS)

    collect_lake_counter = 0

    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    call MAPL_TimerOff(MAPL, "Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")  

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

  ! -------------------------------------------------------------------------------------------
  
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

    ! Pointers to imports (1D tile fields)

    real, pointer :: EVAPOUT(      :) => null()
    real, pointer :: RUNOFF(       :) => null()
    real, pointer :: SHOUT(        :) => null()
    real, pointer :: LWNDSRF(      :) => null()
    real, pointer :: SWNDSRF(      :) => null()
    real, pointer :: HLATN(        :) => null()
    real, pointer :: TST(          :) => null()
    real, pointer :: QST(          :) => null()
    
    ! Pointers to exports ( 1D tile fields)

    real, pointer :: EVAPOUT_enavg(:) => null()
    real, pointer :: RUNOFF_enavg( :) => null()
    real, pointer :: SHOUT_enavg(  :) => null()
    real, pointer :: LWNDSRF_enavg(:) => null()
    real, pointer :: SWNDSRF_enavg(:) => null()
    real, pointer :: HLATN_enavg(  :) => null()
    real, pointer :: TST_enavg(    :) => null()
    real, pointer :: QST_enavg(    :) => null()

    ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_lake_ens_averaging"
    
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Collect_lake")

    ! Get import pointers (1d fields)

    call MAPL_GetPointer(import, EVAPOUT,        'EVAPOUT',        _RC)
    call MAPL_GetPointer(import, RUNOFF ,        'RUNOFF' ,        _RC)
    call MAPL_GetPointer(import, SHOUT  ,        'SHOUT'  ,        _RC)
    call MAPL_GetPointer(import, LWNDSRF,        'LWNDSRF',        _RC)
    call MAPL_GetPointer(import, SWNDSRF,        'SWNDSRF',        _RC)
    call MAPL_GetPointer(import, HLATN  ,        'HLATN'  ,        _RC)
    call MAPL_GetPointer(import, TST    ,        'TST'    ,        _RC)
    call MAPL_GetPointer(import, QST    ,        'QST'    ,        _RC)

    ! Get export pointers (1d fields)

!    call MAPL_GetPointer(import, EVAPOUT_enavg,  'EVAPOUT',        _RC)
!    call MAPL_GetPointer(import, RUNOFF_enavg ,  'RUNOFF' ,        _RC)
!    call MAPL_GetPointer(import, SHOUT_enavg  ,  'SHOUT'  ,        _RC)
!    call MAPL_GetPointer(import, LWNDSRF_enavg,  'LWNDSRF',        _RC)
!    call MAPL_GetPointer(import, SWNDSRF_enavg,  'SWNDSRF',        _RC)
!    call MAPL_GetPointer(import, HLATN_enavg  ,  'HLATN'  ,        _RC)
!    call MAPL_GetPointer(import, TST_enavg    ,  'TST'    ,        _RC)
!    call MAPL_GetPointer(import, QST_enavg    ,  'QST'    ,        _RC)
    call MAPL_GetPointer(export, EVAPOUT_enavg,  'EVAPOUT', _RC)
    call MAPL_GetPointer(export, RUNOFF_enavg ,  'RUNOFF' , _RC)
    call MAPL_GetPointer(export, SHOUT_enavg  ,  'SHOUT'  , _RC)
    call MAPL_GetPointer(export, LWNDSRF_enavg,  'LWNDSRF', _RC)
    call MAPL_GetPointer(export, SWNDSRF_enavg,  'SWNDSRF', _RC)
    call MAPL_GetPointer(export, HLATN_enavg  ,  'HLATN'  , _RC)
    call MAPL_GetPointer(export, TST_enavg    ,  'TST'    , _RC)
    call MAPL_GetPointer(export, QST_enavg    ,  'QST'    , _RC)    
    
    ! On first ensemble member: zero the export accumulators
    
    if (collect_lake_counter == 0) then
       
       if (associated(EVAPOUT_enavg))       EVAPOUT_enavg       = 0.0
       if (associated(RUNOFF_enavg ))       RUNOFF_enavg        = 0.0
       if (associated(SHOUT_enavg  ))       SHOUT_enavg         = 0.0
       if (associated(LWNDSRF_enavg))       LWNDSRF_enavg       = 0.0
       if (associated(SWNDSRF_enavg))       SWNDSRF_enavg       = 0.0
       if (associated(HLATN_enavg  ))       HLATN_enavg         = 0.0
       if (associated(TST_enavg    ))       TST_enavg           = 0.0
       if (associated(QST_enavg    ))       QST_enavg           = 0.0
       
    endif

    ! Accumulate ensemble members (1d fields)
    
    if (associated(EVAPOUT_enavg) .and. associated(EVAPOUT)) EVAPOUT_enavg = EVAPOUT_enavg + EVAPOUT
    if (associated(RUNOFF_enavg)  .and. associated(RUNOFF))  RUNOFF_enavg  = RUNOFF_enavg  + RUNOFF
    if (associated(SHOUT_enavg)   .and. associated(SHOUT))   SHOUT_enavg   = SHOUT_enavg   + SHOUT
    if (associated(LWNDSRF_enavg) .and. associated(LWNDSRF)) LWNDSRF_enavg = LWNDSRF_enavg + LWNDSRF
    if (associated(SWNDSRF_enavg) .and. associated(SWNDSRF)) SWNDSRF_enavg = SWNDSRF_enavg + SWNDSRF
    if (associated(HLATN_enavg)   .and. associated(HLATN))   HLATN_enavg   = HLATN_enavg   + HLATN
    if (associated(TST_enavg)     .and. associated(TST))     TST_enavg     = TST_enavg     + TST
    if (associated(QST_enavg)     .and. associated(QST))     QST_enavg     = QST_enavg     + QST    
    
    collect_lake_counter = collect_lake_counter + 1         ! increment ens member counter 
    
    if (collect_lake_counter == NUM_ENSEMBLE_LAKE) then

       ! reset ens member counter and normalize
       
       collect_lake_counter = 0
       
       ! Divide by NUM_ENSEMBLE_LAKE (1d fields)

       if (associated(EVAPOUT_enavg))       EVAPOUT_enavg       = EVAPOUT_enavg       / NUM_ENSEMBLE_LAKE
       if (associated(RUNOFF_enavg ))       RUNOFF_enavg        = RUNOFF_enavg        / NUM_ENSEMBLE_LAKE
       if (associated(SHOUT_enavg  ))       SHOUT_enavg         = SHOUT_enavg         / NUM_ENSEMBLE_LAKE
       if (associated(LWNDSRF_enavg))       LWNDSRF_enavg       = LWNDSRF_enavg       / NUM_ENSEMBLE_LAKE
       if (associated(SWNDSRF_enavg))       SWNDSRF_enavg       = SWNDSRF_enavg       / NUM_ENSEMBLE_LAKE
       if (associated(HLATN_enavg  ))       HLATN_enavg         = HLATN_enavg         / NUM_ENSEMBLE_LAKE
       if (associated(TST_enavg    ))       TST_enavg           = TST_enavg           / NUM_ENSEMBLE_LAKE
       if (associated(QST_enavg    ))       QST_enavg           = QST_enavg           / NUM_ENSEMBLE_LAKE
       
    endif

    call MAPL_TimerOff(MAPL, "Collect_lake")
    call MAPL_TimerOff(MAPL, "TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

  ! -------------------------------------------------------------------------------------------
  
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
  
end module GEOS_LakeAvgGridCompMod

! ======================= EOF ================================================================================
