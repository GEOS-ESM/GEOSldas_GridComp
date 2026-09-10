#include "MAPL_Generic.h"

module GEOS_MetforceAvgGridCompMod
  use ESMF
  use MAPL_Mod

  use, intrinsic :: ieee_arithmetic
  
  implicit none

  private

  public :: SetServices

  integer            :: NUM_ENSEMBLE
  integer            :: collect_force_counter
  
contains
  
  subroutine SetServices(gc, rc)
    type(ESMF_GridComp), intent(inout) :: gc ! gridded component
    integer, optional                  :: rc ! return code

    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name


    ! Get my name and setup traceback handle
    Iam = 'SetServices'
    call ESMF_GridCompGet(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::" // Iam

    ! Register services for this component
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


!! EXPORT FORCING STATE

   call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "TA",                                                 &
         LONG_NAME  = "perturbed_surface_air_temperature",                      &
         UNITS      = "K",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "QA",                                                 &
         LONG_NAME  = "perturbed_surface_air_specific_humidity",                &
         UNITS      = "kg kg-1",                                                &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PS",                                                 &
         LONG_NAME  = "surface_pressure",                             &
         UNITS      = "Pa",                                                     &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "UU",                                                 &
         LONG_NAME  = "perturbed_surface_wind_speed",                           &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "UWINDLMTILE",                                        &
         LONG_NAME  = "perturbed_levellm_uwind",                                &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "VWINDLMTILE",                                        &
         LONG_NAME  = "perturbed_levellm_vwind",                                &
         UNITS      = "m s-1",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)


    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PCU",                                                &
         LONG_NAME  = "perturbed_liquid_water_convective_precipitation",        &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "PLS",                                                &
         LONG_NAME  = "perturbed_liquid_water_large_scale_precipitation",       &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "RainfSnowf",                                                  &
         LONG_NAME  = "rainf+snowf",                                         &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "SNO",                                                &
         LONG_NAME  = "perturbed_snowfall",                                     &
         UNITS      = "kg m-2 s-1",                                             &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRPAR",                                              &
         LONG_NAME  = "surface_downwelling_PAR_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFPAR",                                              &
         LONG_NAME  = "surface_downwelling_PAR_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

   call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRNIR",                                              &
         LONG_NAME  = "surface_downwelling_nir_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFNIR",                                              &
         LONG_NAME  = "surface_downwelling_nir_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DRUVR",                                              &
         LONG_NAME  = "surface_downwelling_uvr_beam_flux",                      &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DFUVR",                                              &
         LONG_NAME  = "surface_downwelling_uvr_diffuse_flux",                   &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "LWDNSRF",                                            &
         LONG_NAME  = "perturbed_surface_absorbed_longwave_flux",               &
         UNITS      = "W m-2",                                                  &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    call MAPL_AddExportSpec(                                                    &
         gc,                                                                    &
         SHORT_NAME = "DZ",                                                   &
         LONG_NAME  = "reference_height_for_Tair_Qair_Wind",                    &
         UNITS      = "m",                                                      &
         DIMS       = MAPL_DimsTileOnly,                                        &
         VLOCATION  = MAPL_VlocationNone,                                       &
         rc         = status                                                    &
         )
    VERIFY_(status)

    ! Set profiling timers
    call MAPL_TimerAdd(gc, name="Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Collect_force", rc=status)
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

    collect_force_counter    = 0

    call MAPL_GenericInitialize(gc, import, export, clock, rc=status)
    VERIFY_(status)

    call MAPL_TimerOff(MAPL, "Initialize")
    call MAPL_TimerOff(MAPL, "TOTAL")  

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

  subroutine RUN(gc, import, export, clock, rc)
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
    real, pointer :: TApert(:)=>null()
    real, pointer :: QApert(:)=>null()
    real, pointer :: PSpert(:)=>null()
    real, pointer :: UUpert(:)=>null()
    real, pointer :: PCUpert(:)=>null()
    real, pointer :: PLSpert(:)=>null()
    real, pointer :: SNOpert(:)=>null()
    real, pointer :: DRPARpert(:)=>null()
    real, pointer :: DFPARpert(:)=>null()
    real, pointer :: DRNIRpert(:)=>null()
    real, pointer :: DFNIRpert(:)=>null()
    real, pointer :: DRUVRpert(:)=>null()
    real, pointer :: DFUVRpert(:)=>null()
    real, pointer :: LWDNSRFpert(:)=>null()
    real, pointer :: DZpert(:)=>null()

    real, pointer :: TA_enavg(:)=>null()
    real, pointer :: QA_enavg(:)=>null()
    real, pointer :: PS_enavg(:)=>null()
    real, pointer :: UU_enavg(:)=>null()
    real, pointer :: PCU_enavg(:)=>null()
    real, pointer :: PLS_enavg(:)=>null()
    real, pointer :: SNO_enavg(:)=>null()
    real, pointer :: DRPAR_enavg(:)=>null()
    real, pointer :: DFPAR_enavg(:)=>null()
    real, pointer :: DRNIR_enavg(:)=>null()
    real, pointer :: DFNIR_enavg(:)=>null()
    real, pointer :: DRUVR_enavg(:)=>null()
    real, pointer :: DFUVR_enavg(:)=>null()
    real, pointer :: LWDNSRF_enavg(:)=>null()
    real, pointer :: DZ_enavg(:)=>null()
    real, pointer :: RainfSnowf(:)=>null()

   ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_collect_force_ens"

    if(NUM_ENSEMBLE ==1) then
    !   RETURN_(ESMF_SUCCESS)
    endif

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Collect_force")

    call MAPL_GetPointer(import, TApert, 'TApert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, QApert, 'QApert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PSpert, 'PSpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, UUpert, 'UUpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PCUpert, 'PCUpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, PLSpert, 'PLSpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, SNOpert, 'SNOpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRPARpert, 'DRPARpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFPARpert, 'DFPARpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRNIRpert, 'DRNIRpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFNIRpert, 'DFNIRpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DRUVRpert, 'DRUVRpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DFUVRpert, 'DFUVRpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, LWDNSRFpert, 'LWDNSRFpert', rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(import, DZpert, 'DZpert', rc=status)
    VERIFY_(status)

   ! Pointers to exports (allocate memory)
    call MAPL_GetPointer(export, TA_enavg, 'TA', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, QA_enavg, 'QA', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PS_enavg, 'PS', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, UU_enavg, 'UU', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PCU_enavg, 'PCU',alloc=.true.,  rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, PLS_enavg, 'PLS', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, SNO_enavg, 'SNO', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRPAR_enavg, 'DRPAR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFPAR_enavg, 'DFPAR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRNIR_enavg, 'DRNIR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFNIR_enavg, 'DFNIR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DRUVR_enavg, 'DRUVR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DFUVR_enavg, 'DFUVR', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, LWDNSRF_enavg, 'LWDNSRF', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, DZ_enavg, 'DZ', alloc=.true., rc=status)
    VERIFY_(status)
    call MAPL_GetPointer(export, RainfSnowf, 'RainfSnowf', alloc=.true., rc=status)
    VERIFY_(status)

    if(collect_force_counter ==0) then
       if(associated(TA_enavg))  TA_enavg = 0.0
       if(associated(QA_enavg))  QA_enavg = 0.0
       if(associated(PS_enavg))  PS_enavg = 0.0
       if(associated(UU_enavg))  UU_enavg = 0.0
       if(associated(PCU_enavg))  PCU_enavg = 0.0
       if(associated(PLS_enavg))  PLS_enavg = 0.0
       if(associated(SNO_enavg))  SNO_enavg = 0.0
       if(associated(DRPAR_enavg))  DRPAR_enavg = 0.0
       if(associated(DFPAR_enavg))  DFPAR_enavg = 0.0
       if(associated(DRNIR_enavg))  DRNIR_enavg = 0.0
       if(associated(DFNIR_enavg))  DFNIR_enavg = 0.0
       if(associated(DRUVR_enavg))  DRUVR_enavg = 0.0
       if(associated(DFUVR_enavg))  DFUVR_enavg = 0.0
       if(associated(LWDNSRF_enavg))  LWDNSRF_enavg = 0.0
       if(associated(DZ_enavg))  DZ_enavg = 0.0
    endif

    if(associated(TA_enavg))   &
       TA_enavg = TA_enavg + TApert
    if(associated(QA_enavg))   &
       QA_enavg = QA_enavg + QApert
    if(associated(PS_enavg))   &
       PS_enavg = PS_enavg + PSpert
    if(associated(UU_enavg))   &
       UU_enavg = UU_enavg + UUpert
    if(associated(PCU_enavg))   &
       PCU_enavg = PCU_enavg + PCUpert
    if(associated(PLS_enavg))   &
       PLS_enavg = PLS_enavg + PLSpert
    if(associated(SNO_enavg))   &
       SNO_enavg = SNO_enavg + SNOpert
    if(associated(DRPAR_enavg))   &
       DRPAR_enavg = DRPAR_enavg + DRPARpert
    if(associated(DFPAR_enavg))   &
       DFPAR_enavg = DFPAR_enavg + DFPARpert
    if(associated(DRNIR_enavg))   &
       DRNIR_enavg = DRNIR_enavg + DRNIRpert
    if(associated(DFNIR_enavg))   &
       DFNIR_enavg = DFNIR_enavg + DFNIRpert
    if(associated(DRUVR_enavg))   &
       DRUVR_enavg = DRUVR_enavg + DRUVRpert
    if(associated(DFUVR_enavg))   &
       DFUVR_enavg = DFUVR_enavg + DFUVRpert
    if(associated(LWDNSRF_enavg))   &
       LWDNSRF_enavg = LWDNSRF_enavg + LWDNSRFpert
    if(associated(DZ_enavg))   &
       DZ_enavg = DZ_enavg + DZpert

    collect_force_counter = collect_force_counter + 1
    if(collect_force_counter == NUM_ENSEMBLE) then
      collect_force_counter = 0

      if(associated(TA_enavg))  TA_enavg =TA_enavg /NUM_ENSEMBLE
      if(associated(QA_enavg))  QA_enavg =QA_enavg /NUM_ENSEMBLE
      if(associated(PS_enavg))  PS_enavg =PS_enavg /NUM_ENSEMBLE
      if(associated(UU_enavg))  UU_enavg =UU_enavg /NUM_ENSEMBLE
      if(associated(PCU_enavg))  PCU_enavg =PCU_enavg /NUM_ENSEMBLE
      if(associated(PLS_enavg))  PLS_enavg =PLS_enavg /NUM_ENSEMBLE
      if(associated(SNO_enavg))  SNO_enavg =SNO_enavg /NUM_ENSEMBLE
      if(associated(DRPAR_enavg))  DRPAR_enavg =DRPAR_enavg /NUM_ENSEMBLE
      if(associated(DFPAR_enavg))  DFPAR_enavg =DFPAR_enavg /NUM_ENSEMBLE
      if(associated(DRNIR_enavg))  DRNIR_enavg =DRNIR_enavg /NUM_ENSEMBLE
      if(associated(DFNIR_enavg))  DFNIR_enavg =DFNIR_enavg /NUM_ENSEMBLE
      if(associated(DRUVR_enavg))  DRUVR_enavg =DRUVR_enavg /NUM_ENSEMBLE
      if(associated(DFUVR_enavg))  DFUVR_enavg =DFUVR_enavg /NUM_ENSEMBLE
      if(associated(LWDNSRF_enavg))  LWDNSRF_enavg =LWDNSRF_enavg /NUM_ENSEMBLE
      if(associated(DZ_enavg))   DZ_enavg      =DZ_enavg /NUM_ENSEMBLE
      if(associated(RainfSnowf)) RainfSnowf    =PLS_enavg + PCU_enavg + SNO_enavg

    endif
    ! Turn timers off
    call MAPL_TimerOff(MAPL, "Collect_force")
    call MAPL_TimerOff(MAPL, "TOTAL")
       
   RETURN_(ESMF_SUCCESS)

  end subroutine RUN

  subroutine Finalize(gc, import, export, clock, rc)

    ! !ARGUMENTS:

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

end module GEOS_MetforceAvgGridCompMod
