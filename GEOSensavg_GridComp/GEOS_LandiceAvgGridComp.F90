#include "MAPL_Generic.h"

module GEOS_LandiceAvgGridCompMod
  use ESMF
  use MAPL_Mod
  use GEOS_LandiceGridCompMod, only: NUM_SNOW_LAYERS,  NUM_ICE_LAYERS
  implicit none

  private

  public   :: SetServices

  integer  :: NUM_ENSEMBLE
  integer  :: collect_landice_counter

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

!  export landice

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_emissivity'                              ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'EMIS'                      ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_reflectivity_for_visible_beam'           ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'ALBVR'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_reflectivity_for_visible_diffuse'        ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'ALBVF'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_reflectivity_for_near_infrared_beam'     ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'ALBNR'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_reflectivity_for_near_infrared_diffuse'  ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'ALBNF'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_temperature'                             ,&
         UNITS              = 'K'                         ,&
         SHORT_NAME         = 'TST'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'land_surface_skin_temperature'                   ,&
         UNITS              = 'K'                         ,&
         SHORT_NAME         = 'LST'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_specific_humidity'                       ,&
         UNITS              = 'kg kg-1'                   ,&
         SHORT_NAME         = 'QST'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'turbulence_surface_skin_temperature'             ,&
         UNITS              = 'K'                         ,&
         SHORT_NAME         = 'TH'                        ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'turbulence_surface_specific_humidity'            ,&
         UNITS              = 'kg kg-1'                   ,&
         SHORT_NAME         = 'QH'                        ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'change_of_surface_skin_temperature'              ,&
         UNITS              = 'K'                         ,&
         SHORT_NAME         = 'DELTS'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'change_of_surface_specific_humidity'             ,&
         UNITS              = 'kg kg-1'                   ,&
         SHORT_NAME         = 'DELQS'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_heat_exchange_coefficient'               ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'CHT'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_momentum_exchange_coefficient'           ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'CMT'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_moisture_exchange_coefficient'           ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'CQT'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'neutral_drag_coefficient'                        ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'CNT'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_bulk_richardson_number'                  ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'RIT'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'net_ice_accumulation_rate'                       ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'ACCUM'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_ice_evaporation_energy_flux_over_glaciated_surface',&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'EVPICE_GL'                 ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'sublimation'                                     ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'SUBLIM'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_mass_over_glaciated_surface'                ,&
         UNITS              = 'kg m-2'                    ,&
         SHORT_NAME         = 'SNOMAS_GL'                 ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_mass_over_glaciated_surface'                ,&
         UNITS              = 'kg m-2'                    ,&
         SHORT_NAME         = 'SNOWMASS'                  ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_depth_over_glaciated_surface'               ,&
         UNITS              = 'm'                         ,&
         SHORT_NAME         = 'SNOWDP_GL'                 ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'fractional_snow_covered_area_of_glaciated_surface',&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'ASNOW_GL'                  ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_density'                              ,&
         UNITS              = 'kg m-3'                    ,&
         SHORT_NAME         = 'RHOSNOW'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_temperature'                          ,&
         UNITS              = 'deg C'                     ,&
         SHORT_NAME         = 'TSNOW'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'aggregated_ice_layer_temperature'                ,&
         UNITS              = 'deg C'                     ,&
         SHORT_NAME         = 'TICE0'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_ICE_LAYERS/)          ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_water_content'                        ,&
         UNITS              = 'kg m-2'                    ,&
         SHORT_NAME         = 'WSNOW'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_thickness'                            ,&
         UNITS              = 'm'                         ,&
         SHORT_NAME         = 'ZSNOW'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_density_change_due_to_densification'  ,&
         UNITS              = 'kg m-3'                    ,&
         SHORT_NAME         = 'DRHOS0'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_mass_residual_due_to_densification'   ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'WESNEX'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'total_snow_mass_residual_due_to_densification'   ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'WESNEXT'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'top_snow_layer_mass_change_due_to_sublimation_and_condensation',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'WESNSC'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'top_snow_layer_thickness_change_due_to_sub_con'  ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'SNDZSC'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'top_snow_layer_mass_change_due_to_precip'        ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'WESNPREC'                  ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'top_snow_layer_thickness_change_due_to_precip'   ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'SNDZPREC'                  ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'top_snow_layer_thickness_change_due_to_percolation',&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'SNDZ1PERC'                 ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_mass_change_due_to_percolation'       ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'WESNPERC'                  ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_mass_change_due_to_densification'     ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'WESNDENS'                  ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_layer_mass_change_due_to_repartition'       ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'WESNREPAR'                 ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS/)         ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'frozen_runoff_due_to_fixed_max_depth'            ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'WESNBOT'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'contribution_to_surface_mass_balance_from_rain_frozen_onto_bare_ice',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RAINRFZ'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_melt_flux'                                  ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'SMELT'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'ice_melt_flux'                                   ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'IMELT'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snow_broadband_reflectivity'                     ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'SNOWALB'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'aggregated_snow_ice_broadband_reflectivity'      ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'SNICEALB'                  ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'melt_water_production'                           ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'MELTWTR'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'snowpack_meltwater_content'                      ,&
         UNITS              = 'kg m-2'                    ,&
         SHORT_NAME         = 'MELTWTRCONT'               ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'liquid_water_content_in_top_x_m'                ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'LWC'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'runoff_total_flux'                               ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RUNOFF'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'gustiness'                                       ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'GUST'                      ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_ventilation_velocity'                    ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'VENT'                      ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_roughness'                               ,&
         UNITS              = 'm'                         ,&
         SHORT_NAME         = 'Z0'                        ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_roughness_for_heat'                      ,&
         UNITS              = 'm'                         ,&
         SHORT_NAME         = 'Z0H'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'temperature_2m_wind_from_MO_sfc'                ,&
         UNITS              = 'K'                         ,&
         SHORT_NAME         = 'MOT2M'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'humidity_2m_wind_from_MO_sfc'                   ,&
         UNITS              = 'kg kg-1'                   ,&
         SHORT_NAME         = 'MOQ2M'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'zonal_2m_wind_from_MO_sfc'                      ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'MOU2M'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'meridional_2m_wind_from_MO_sfc'                 ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'MOV2M'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'temperature_10m_wind_from_MO_sfc'               ,&
         UNITS              = 'K'                         ,&
         SHORT_NAME         = 'MOT10M'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'humidity_10m_wind_from_MO_sfc'                  ,&
         UNITS              = 'kg kg-1'                   ,&
         SHORT_NAME         = 'MOQ10M'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'zonal_10m_wind_from_MO_sfc'                     ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'MOU10M'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'meridional_10m_wind_from_MO_sfc'                ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'MOV10M'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'zonal_50m_wind_from_MO_sfc'                     ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'MOU50M'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'meridional_50m_wind_from_MO_sfc'                ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'MOV50M'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'evaporation'                                     ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'EVAPOUT'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'upward_sensible_heat_flux'                       ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'SHOUT'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_emitted_longwave_flux'                   ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'HLWUP'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_net_downward_longwave_flux'              ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'LWNDSRF'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'surface_net_downward_shortwave_flux'             ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'SWNDSRF'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'total_latent_energy_flux'                        ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'HLATN'                     ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'downward_heat_flux_in_ice'                       ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'DNICFLX'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'Ground_heating_snow'                             ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'GHSNOW'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'glacier_ice_heating_flux'                        ,&
         UNITS              = 'W m-2'                     ,&
         SHORT_NAME         = 'GHTSKIN'                   ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'vegetation_type'                                 ,&
         UNITS              = '1'                         ,&
         SHORT_NAME         = 'ITY'                       ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_1',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTDU001'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_2',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTDU002'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_3',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTDU003'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_4',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTDU004'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_5',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTDU005'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_1',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTBC001'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_2',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTBC002'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_1',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTOC001'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_2',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RMELTOC002'                ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    ! Set profiling timers
    call MAPL_TimerAdd(gc, name="Initialize", rc=status)
    VERIFY_(status)
    call MAPL_TimerAdd(gc, name="Collect_landice", rc=status)
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

    collect_landice_counter = 0

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

    ! Pointers to imports (1D tile fields)
    real, pointer :: EMIS(:)=>null()
    real, pointer :: ALBVR(:)=>null()
    real, pointer :: ALBVF(:)=>null()
    real, pointer :: ALBNR(:)=>null()
    real, pointer :: ALBNF(:)=>null()
    real, pointer :: TST(:)=>null()
    real, pointer :: LST(:)=>null()
    real, pointer :: QST(:)=>null()
    real, pointer :: TH(:)=>null()
    real, pointer :: QH(:)=>null()
    real, pointer :: DELTS(:)=>null()
    real, pointer :: DELQS(:)=>null()
    real, pointer :: CHT(:)=>null()
    real, pointer :: CMT(:)=>null()
    real, pointer :: CQT(:)=>null()
    real, pointer :: CNT(:)=>null()
    real, pointer :: RIT(:)=>null()
    real, pointer :: ACCUM(:)=>null()
    real, pointer :: EVPICE_GL(:)=>null()
    real, pointer :: SUBLIM(:)=>null()
    real, pointer :: SNOMAS_GL(:)=>null()
    real, pointer :: SNOWMASS(:)=>null()
    real, pointer :: SNOWDP_GL(:)=>null()
    real, pointer :: ASNOW_GL(:)=>null()
    real, pointer :: WESNEXT(:)=>null()
    real, pointer :: WESNSC(:)=>null()
    real, pointer :: SNDZSC(:)=>null()
    real, pointer :: WESNPREC(:)=>null()
    real, pointer :: SNDZPREC(:)=>null()
    real, pointer :: SNDZ1PERC(:)=>null()
    real, pointer :: WESNBOT(:)=>null()
    real, pointer :: RAINRFZ(:)=>null()
    real, pointer :: SMELT(:)=>null()
    real, pointer :: IMELT(:)=>null()
    real, pointer :: SNOWALB(:)=>null()
    real, pointer :: SNICEALB(:)=>null()
    real, pointer :: MELTWTR(:)=>null()
    real, pointer :: MELTWTRCONT(:)=>null()
    real, pointer :: LWC(:)=>null()
    real, pointer :: RUNOFF(:)=>null()
    real, pointer :: GUST(:)=>null()
    real, pointer :: VENT(:)=>null()
    real, pointer :: Z0(:)=>null()
    real, pointer :: Z0H(:)=>null()
    real, pointer :: MOT2M(:)=>null()
    real, pointer :: MOQ2M(:)=>null()
    real, pointer :: MOU2M(:)=>null()
    real, pointer :: MOV2M(:)=>null()
    real, pointer :: MOT10M(:)=>null()
    real, pointer :: MOQ10M(:)=>null()
    real, pointer :: MOU10M(:)=>null()
    real, pointer :: MOV10M(:)=>null()
    real, pointer :: MOU50M(:)=>null()
    real, pointer :: MOV50M(:)=>null()
    real, pointer :: EVAPOUT(:)=>null()
    real, pointer :: SHOUT(:)=>null()
    real, pointer :: HLWUP(:)=>null()
    real, pointer :: LWNDSRF(:)=>null()
    real, pointer :: SWNDSRF(:)=>null()
    real, pointer :: HLATN(:)=>null()
    real, pointer :: DNICFLX(:)=>null()
    real, pointer :: GHSNOW(:)=>null()
    real, pointer :: GHTSKIN(:)=>null()
    real, pointer :: ITY(:)=>null()
    real, pointer :: RMELTDU001(:)=>null()
    real, pointer :: RMELTDU002(:)=>null()
    real, pointer :: RMELTDU003(:)=>null()
    real, pointer :: RMELTDU004(:)=>null()
    real, pointer :: RMELTDU005(:)=>null()
    real, pointer :: RMELTBC001(:)=>null()
    real, pointer :: RMELTBC002(:)=>null()
    real, pointer :: RMELTOC001(:)=>null()
    real, pointer :: RMELTOC002(:)=>null()

    ! Pointers to imports (multi-dimensional fields with UNGRIDDED_DIMS)
    real, pointer :: RHOSNOW(:,:)=>null()
    real, pointer :: TSNOW(:,:)=>null()
    real, pointer :: TICE0(:,:)=>null()
    real, pointer :: WSNOW(:,:)=>null()
    real, pointer :: ZSNOW(:,:)=>null()
    real, pointer :: DRHOS0(:,:)=>null()
    real, pointer :: WESNEX(:,:)=>null()
    real, pointer :: WESNPERC(:,:)=>null()
    real, pointer :: WESNDENS(:,:)=>null()
    real, pointer :: WESNREPAR(:,:)=>null()

    ! Pointers to exports ( 1D tile fields)
    real, pointer :: EMIS_enavg(:)=>null()
    real, pointer :: ALBVR_enavg(:)=>null()
    real, pointer :: ALBVF_enavg(:)=>null()
    real, pointer :: ALBNR_enavg(:)=>null()
    real, pointer :: ALBNF_enavg(:)=>null()
    real, pointer :: TST_enavg(:)=>null()
    real, pointer :: LST_enavg(:)=>null()
    real, pointer :: QST_enavg(:)=>null()
    real, pointer :: TH_enavg(:)=>null()
    real, pointer :: QH_enavg(:)=>null()
    real, pointer :: DELTS_enavg(:)=>null()
    real, pointer :: DELQS_enavg(:)=>null()
    real, pointer :: CHT_enavg(:)=>null()
    real, pointer :: CMT_enavg(:)=>null()
    real, pointer :: CQT_enavg(:)=>null()
    real, pointer :: CNT_enavg(:)=>null()
    real, pointer :: RIT_enavg(:)=>null()
    real, pointer :: ACCUM_enavg(:)=>null()
    real, pointer :: EVPICE_GL_enavg(:)=>null()
    real, pointer :: SUBLIM_enavg(:)=>null()
    real, pointer :: SNOMAS_GL_enavg(:)=>null()
    real, pointer :: SNOWMASS_enavg(:)=>null()
    real, pointer :: SNOWDP_GL_enavg(:)=>null()
    real, pointer :: ASNOW_GL_enavg(:)=>null()
    real, pointer :: WESNEXT_enavg(:)=>null()
    real, pointer :: WESNSC_enavg(:)=>null()
    real, pointer :: SNDZSC_enavg(:)=>null()
    real, pointer :: WESNPREC_enavg(:)=>null()
    real, pointer :: SNDZPREC_enavg(:)=>null()
    real, pointer :: SNDZ1PERC_enavg(:)=>null()
    real, pointer :: WESNBOT_enavg(:)=>null()
    real, pointer :: RAINRFZ_enavg(:)=>null()
    real, pointer :: SMELT_enavg(:)=>null()
    real, pointer :: IMELT_enavg(:)=>null()
    real, pointer :: SNOWALB_enavg(:)=>null()
    real, pointer :: SNICEALB_enavg(:)=>null()
    real, pointer :: MELTWTR_enavg(:)=>null()
    real, pointer :: MELTWTRCONT_enavg(:)=>null()
    real, pointer :: LWC_enavg(:)=>null()
    real, pointer :: RUNOFF_enavg(:)=>null()
    real, pointer :: GUST_enavg(:)=>null()
    real, pointer :: VENT_enavg(:)=>null()
    real, pointer :: Z0_enavg(:)=>null()
    real, pointer :: Z0H_enavg(:)=>null()
    real, pointer :: MOT2M_enavg(:)=>null()
    real, pointer :: MOQ2M_enavg(:)=>null()
    real, pointer :: MOU2M_enavg(:)=>null()
    real, pointer :: MOV2M_enavg(:)=>null()
    real, pointer :: MOT10M_enavg(:)=>null()
    real, pointer :: MOQ10M_enavg(:)=>null()
    real, pointer :: MOU10M_enavg(:)=>null()
    real, pointer :: MOV10M_enavg(:)=>null()
    real, pointer :: MOU50M_enavg(:)=>null()
    real, pointer :: MOV50M_enavg(:)=>null()
    real, pointer :: EVAPOUT_enavg(:)=>null()
    real, pointer :: SHOUT_enavg(:)=>null()
    real, pointer :: HLWUP_enavg(:)=>null()
    real, pointer :: LWNDSRF_enavg(:)=>null()
    real, pointer :: SWNDSRF_enavg(:)=>null()
    real, pointer :: HLATN_enavg(:)=>null()
    real, pointer :: DNICFLX_enavg(:)=>null()
    real, pointer :: GHSNOW_enavg(:)=>null()
    real, pointer :: GHTSKIN_enavg(:)=>null()
    real, pointer :: ITY_enavg(:)=>null()
    real, pointer :: RMELTDU001_enavg(:)=>null()
    real, pointer :: RMELTDU002_enavg(:)=>null()
    real, pointer :: RMELTDU003_enavg(:)=>null()
    real, pointer :: RMELTDU004_enavg(:)=>null()
    real, pointer :: RMELTDU005_enavg(:)=>null()
    real, pointer :: RMELTBC001_enavg(:)=>null()
    real, pointer :: RMELTBC002_enavg(:)=>null()
    real, pointer :: RMELTOC001_enavg(:)=>null()
    real, pointer :: RMELTOC002_enavg(:)=>null()

    ! Pointers to exports (multi-dimensional fields with UNGRIDDED_DIMS)
    real, pointer :: RHOSNOW_enavg(:,:)=>null()
    real, pointer :: TSNOW_enavg(:,:)=>null()
    real, pointer :: TICE0_enavg(:,:)=>null()
    real, pointer :: WSNOW_enavg(:,:)=>null()
    real, pointer :: ZSNOW_enavg(:,:)=>null()
    real, pointer :: DRHOS0_enavg(:,:)=>null()
    real, pointer :: WESNEX_enavg(:,:)=>null()
    real, pointer :: WESNPERC_enavg(:,:)=>null()
    real, pointer :: WESNDENS_enavg(:,:)=>null()
    real, pointer :: WESNREPAR_enavg(:,:)=>null()

     ! Get my name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Run_landice_ens_averaging"

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    ! Turn timers on
    call MAPL_TimerOn(MAPL, "TOTAL")
    call MAPL_TimerOn(MAPL, "Collect_landice")

    ! Get import pointers (1d fields)
    call MAPL_GetPointer(import, EMIS,        'EMIS',        _RC)
    call MAPL_GetPointer(import, ALBVR,       'ALBVR',       _RC)
    call MAPL_GetPointer(import, ALBVF,       'ALBVF',       _RC)
    call MAPL_GetPointer(import, ALBNR,       'ALBNR',       _RC)
    call MAPL_GetPointer(import, ALBNF,       'ALBNF',       _RC)
    call MAPL_GetPointer(import, TST,         'TST',         _RC)
    call MAPL_GetPointer(import, LST,         'LST',         _RC)
    call MAPL_GetPointer(import, QST,         'QST',         _RC)
    call MAPL_GetPointer(import, TH,          'TH',          _RC)
    call MAPL_GetPointer(import, QH,          'QH',          _RC)
    call MAPL_GetPointer(import, DELTS,       'DELTS',       _RC)
    call MAPL_GetPointer(import, DELQS,       'DELQS',       _RC)
    call MAPL_GetPointer(import, CHT,         'CHT',         _RC)
    call MAPL_GetPointer(import, CMT,         'CMT',         _RC)
    call MAPL_GetPointer(import, CQT,         'CQT',         _RC)
    call MAPL_GetPointer(import, CNT,         'CNT',         _RC)
    call MAPL_GetPointer(import, RIT,         'RIT',         _RC)
    call MAPL_GetPointer(import, ACCUM,       'ACCUM',       _RC)
    call MAPL_GetPointer(import, EVPICE_GL,   'EVPICE_GL',   _RC)
    call MAPL_GetPointer(import, SUBLIM,      'SUBLIM',      _RC)
    call MAPL_GetPointer(import, SNOMAS_GL,   'SNOMAS_GL',   _RC)
    call MAPL_GetPointer(import, SNOWMASS,    'SNOWMASS',    _RC)
    call MAPL_GetPointer(import, SNOWDP_GL,   'SNOWDP_GL',   _RC)
    call MAPL_GetPointer(import, ASNOW_GL,    'ASNOW_GL',    _RC)
    call MAPL_GetPointer(import, WESNEXT,     'WESNEXT',     _RC)
    call MAPL_GetPointer(import, WESNSC,      'WESNSC',      _RC)
    call MAPL_GetPointer(import, SNDZSC,      'SNDZSC',      _RC)
    call MAPL_GetPointer(import, WESNPREC,    'WESNPREC',    _RC)
    call MAPL_GetPointer(import, SNDZPREC,    'SNDZPREC',    _RC)
    call MAPL_GetPointer(import, SNDZ1PERC,   'SNDZ1PERC',   _RC)
    call MAPL_GetPointer(import, WESNBOT,     'WESNBOT',     _RC)
    call MAPL_GetPointer(import, RAINRFZ,     'RAINRFZ',     _RC)
    call MAPL_GetPointer(import, SMELT,       'SMELT',       _RC)
    call MAPL_GetPointer(import, IMELT,       'IMELT',       _RC)
    call MAPL_GetPointer(import, SNOWALB,     'SNOWALB',     _RC)
    call MAPL_GetPointer(import, SNICEALB,    'SNICEALB',    _RC)
    call MAPL_GetPointer(import, MELTWTR,     'MELTWTR',     _RC)
    call MAPL_GetPointer(import, MELTWTRCONT, 'MELTWTRCONT', _RC)
    call MAPL_GetPointer(import, LWC,         'LWC',         _RC)
    call MAPL_GetPointer(import, RUNOFF,      'RUNOFF',      _RC)
    call MAPL_GetPointer(import, GUST,        'GUST',        _RC)
    call MAPL_GetPointer(import, VENT,        'VENT',        _RC)
    call MAPL_GetPointer(import, Z0,          'Z0',          _RC)
    call MAPL_GetPointer(import, Z0H,         'Z0H',         _RC)
    call MAPL_GetPointer(import, MOT2M,       'MOT2M',       _RC)
    call MAPL_GetPointer(import, MOQ2M,       'MOQ2M',       _RC)
    call MAPL_GetPointer(import, MOU2M,       'MOU2M',       _RC)
    call MAPL_GetPointer(import, MOV2M,       'MOV2M',       _RC)
    call MAPL_GetPointer(import, MOT10M,      'MOT10M',      _RC)
    call MAPL_GetPointer(import, MOQ10M,      'MOQ10M',      _RC)
    call MAPL_GetPointer(import, MOU10M,      'MOU10M',      _RC)
    call MAPL_GetPointer(import, MOV10M,      'MOV10M',      _RC)
    call MAPL_GetPointer(import, MOU50M,      'MOU50M',      _RC)
    call MAPL_GetPointer(import, MOV50M,      'MOV50M',      _RC)
    call MAPL_GetPointer(import, EVAPOUT,     'EVAPOUT',     _RC)
    call MAPL_GetPointer(import, SHOUT,       'SHOUT',       _RC)
    call MAPL_GetPointer(import, HLWUP,       'HLWUP',       _RC)
    call MAPL_GetPointer(import, LWNDSRF,     'LWNDSRF',     _RC)
    call MAPL_GetPointer(import, SWNDSRF,     'SWNDSRF',     _RC)
    call MAPL_GetPointer(import, HLATN,       'HLATN',       _RC)
    call MAPL_GetPointer(import, DNICFLX,     'DNICFLX',     _RC)
    call MAPL_GetPointer(import, GHSNOW,      'GHSNOW',      _RC)
    call MAPL_GetPointer(import, GHTSKIN,     'GHTSKIN',     _RC)
    call MAPL_GetPointer(import, ITY,         'ITY',         _RC)
    call MAPL_GetPointer(import, RMELTDU001,  'RMELTDU001',  _RC)
    call MAPL_GetPointer(import, RMELTDU002,  'RMELTDU002',  _RC)
    call MAPL_GetPointer(import, RMELTDU003,  'RMELTDU003',  _RC)
    call MAPL_GetPointer(import, RMELTDU004,  'RMELTDU004',  _RC)
    call MAPL_GetPointer(import, RMELTDU005,  'RMELTDU005',  _RC)
    call MAPL_GetPointer(import, RMELTBC001,  'RMELTBC001',  _RC)
    call MAPL_GetPointer(import, RMELTBC002,  'RMELTBC002',  _RC)
    call MAPL_GetPointer(import, RMELTOC001,  'RMELTOC001',  _RC)
    call MAPL_GetPointer(import, RMELTOC002,  'RMELTOC002',  _RC)

    ! Get import pointers (multi-dimensional fields)
    call MAPL_GetPointer(import, RHOSNOW,     'RHOSNOW',     _RC)
    call MAPL_GetPointer(import, TSNOW,       'TSNOW',       _RC)
    call MAPL_GetPointer(import, TICE0,       'TICE0',       _RC)
    call MAPL_GetPointer(import, WSNOW,       'WSNOW',       _RC)
    call MAPL_GetPointer(import, ZSNOW,       'ZSNOW',       _RC)
    call MAPL_GetPointer(import, DRHOS0,      'DRHOS0',      _RC)
    call MAPL_GetPointer(import, WESNEX,      'WESNEX',      _RC)
    call MAPL_GetPointer(import, WESNPERC,    'WESNPERC',    _RC)
    call MAPL_GetPointer(import, WESNDENS,    'WESNDENS',    _RC)
    call MAPL_GetPointer(import, WESNREPAR,   'WESNREPAR',   _RC)

    ! Get export pointers (1d fields)
    call MAPL_GetPointer(export, EMIS_enavg,        'EMIS',        _RC)
    call MAPL_GetPointer(export, ALBVR_enavg,       'ALBVR',       _RC)
    call MAPL_GetPointer(export, ALBVF_enavg,       'ALBVF',       _RC)
    call MAPL_GetPointer(export, ALBNR_enavg,       'ALBNR',       _RC)
    call MAPL_GetPointer(export, ALBNF_enavg,       'ALBNF',       _RC)
    call MAPL_GetPointer(export, TST_enavg,         'TST',         _RC)
    call MAPL_GetPointer(export, LST_enavg,         'LST',         _RC)
    call MAPL_GetPointer(export, QST_enavg,         'QST',         _RC)
    call MAPL_GetPointer(export, TH_enavg,          'TH',          _RC)
    call MAPL_GetPointer(export, QH_enavg,          'QH',          _RC)
    call MAPL_GetPointer(export, DELTS_enavg,       'DELTS',       _RC)
    call MAPL_GetPointer(export, DELQS_enavg,       'DELQS',       _RC)
    call MAPL_GetPointer(export, CHT_enavg,         'CHT',         _RC)
    call MAPL_GetPointer(export, CMT_enavg,         'CMT',         _RC)
    call MAPL_GetPointer(export, CQT_enavg,         'CQT',         _RC)
    call MAPL_GetPointer(export, CNT_enavg,         'CNT',         _RC)
    call MAPL_GetPointer(export, RIT_enavg,         'RIT',         _RC)
    call MAPL_GetPointer(export, ACCUM_enavg,       'ACCUM',       _RC)
    call MAPL_GetPointer(export, EVPICE_GL_enavg,   'EVPICE_GL',   _RC)
    call MAPL_GetPointer(export, SUBLIM_enavg,      'SUBLIM',      _RC)
    call MAPL_GetPointer(export, SNOMAS_GL_enavg,   'SNOMAS_GL',   _RC)
    call MAPL_GetPointer(export, SNOWMASS_enavg,    'SNOWMASS',    _RC)
    call MAPL_GetPointer(export, SNOWDP_GL_enavg,   'SNOWDP_GL',   _RC)
    call MAPL_GetPointer(export, ASNOW_GL_enavg,    'ASNOW_GL',    _RC)
    call MAPL_GetPointer(export, WESNEXT_enavg,     'WESNEXT',     _RC)
    call MAPL_GetPointer(export, WESNSC_enavg,      'WESNSC',      _RC)
    call MAPL_GetPointer(export, SNDZSC_enavg,      'SNDZSC',      _RC)
    call MAPL_GetPointer(export, WESNPREC_enavg,    'WESNPREC',    _RC)
    call MAPL_GetPointer(export, SNDZPREC_enavg,    'SNDZPREC',    _RC)
    call MAPL_GetPointer(export, SNDZ1PERC_enavg,   'SNDZ1PERC',   _RC)
    call MAPL_GetPointer(export, WESNBOT_enavg,     'WESNBOT',     _RC)
    call MAPL_GetPointer(export, RAINRFZ_enavg,     'RAINRFZ',     _RC)
    call MAPL_GetPointer(export, SMELT_enavg,       'SMELT',       _RC)
    call MAPL_GetPointer(export, IMELT_enavg,       'IMELT',       _RC)
    call MAPL_GetPointer(export, SNOWALB_enavg,     'SNOWALB',     _RC)
    call MAPL_GetPointer(export, SNICEALB_enavg,    'SNICEALB',    _RC)
    call MAPL_GetPointer(export, MELTWTR_enavg,     'MELTWTR',     _RC)
    call MAPL_GetPointer(export, MELTWTRCONT_enavg, 'MELTWTRCONT', _RC)
    call MAPL_GetPointer(export, LWC_enavg,         'LWC',         _RC)
    call MAPL_GetPointer(export, RUNOFF_enavg,      'RUNOFF',      _RC)
    call MAPL_GetPointer(export, GUST_enavg,        'GUST',        _RC)
    call MAPL_GetPointer(export, VENT_enavg,        'VENT',        _RC)
    call MAPL_GetPointer(export, Z0_enavg,          'Z0',          _RC)
    call MAPL_GetPointer(export, Z0H_enavg,         'Z0H',         _RC)
    call MAPL_GetPointer(export, MOT2M_enavg,       'MOT2M',       _RC)
    call MAPL_GetPointer(export, MOQ2M_enavg,       'MOQ2M',       _RC)
    call MAPL_GetPointer(export, MOU2M_enavg,       'MOU2M',       _RC)
    call MAPL_GetPointer(export, MOV2M_enavg,       'MOV2M',       _RC)
    call MAPL_GetPointer(export, MOT10M_enavg,      'MOT10M',      _RC)
    call MAPL_GetPointer(export, MOQ10M_enavg,      'MOQ10M',      _RC)
    call MAPL_GetPointer(export, MOU10M_enavg,      'MOU10M',      _RC)
    call MAPL_GetPointer(export, MOV10M_enavg,      'MOV10M',      _RC)
    call MAPL_GetPointer(export, MOU50M_enavg,      'MOU50M',      _RC)
    call MAPL_GetPointer(export, MOV50M_enavg,      'MOV50M',      _RC)
    call MAPL_GetPointer(export, EVAPOUT_enavg,     'EVAPOUT',     _RC)
    call MAPL_GetPointer(export, SHOUT_enavg,       'SHOUT',       _RC)
    call MAPL_GetPointer(export, HLWUP_enavg,       'HLWUP',       _RC)
    call MAPL_GetPointer(export, LWNDSRF_enavg,     'LWNDSRF',     _RC)
    call MAPL_GetPointer(export, SWNDSRF_enavg,     'SWNDSRF',     _RC)
    call MAPL_GetPointer(export, HLATN_enavg,       'HLATN',       _RC)
    call MAPL_GetPointer(export, DNICFLX_enavg,     'DNICFLX',     _RC)
    call MAPL_GetPointer(export, GHSNOW_enavg,      'GHSNOW',      _RC)
    call MAPL_GetPointer(export, GHTSKIN_enavg,     'GHTSKIN',     _RC)
    call MAPL_GetPointer(export, ITY_enavg,         'ITY',         _RC)
    call MAPL_GetPointer(export, RMELTDU001_enavg,  'RMELTDU001',  _RC)
    call MAPL_GetPointer(export, RMELTDU002_enavg,  'RMELTDU002',  _RC)
    call MAPL_GetPointer(export, RMELTDU003_enavg,  'RMELTDU003',  _RC)
    call MAPL_GetPointer(export, RMELTDU004_enavg,  'RMELTDU004',  _RC)
    call MAPL_GetPointer(export, RMELTDU005_enavg,  'RMELTDU005',  _RC)
    call MAPL_GetPointer(export, RMELTBC001_enavg,  'RMELTBC001',  _RC)
    call MAPL_GetPointer(export, RMELTBC002_enavg,  'RMELTBC002',  _RC)
    call MAPL_GetPointer(export, RMELTOC001_enavg,  'RMELTOC001',  _RC)
    call MAPL_GetPointer(export, RMELTOC002_enavg,  'RMELTOC002',  _RC)

    ! Get export pointers (multi-dimensional fields)
    call MAPL_GetPointer(export, RHOSNOW_enavg,     'RHOSNOW',     _RC)
    call MAPL_GetPointer(export, TSNOW_enavg,       'TSNOW',       _RC)
    call MAPL_GetPointer(export, TICE0_enavg,       'TICE0',       _RC)
    call MAPL_GetPointer(export, WSNOW_enavg,       'WSNOW',       _RC)
    call MAPL_GetPointer(export, ZSNOW_enavg,       'ZSNOW',       _RC)
    call MAPL_GetPointer(export, DRHOS0_enavg,      'DRHOS0',      _RC)
    call MAPL_GetPointer(export, WESNEX_enavg,      'WESNEX',      _RC)
    call MAPL_GetPointer(export, WESNPERC_enavg,    'WESNPERC',    _RC)
    call MAPL_GetPointer(export, WESNDENS_enavg,    'WESNDENS',    _RC)
    call MAPL_GetPointer(export, WESNREPAR_enavg,   'WESNREPAR',   _RC)

    ! On first ensemble member: zero the export accumulators
    if (collect_landice_counter == 0) then
       if (associated(EMIS_enavg))        EMIS_enavg        = 0.0
       if (associated(ALBVR_enavg))       ALBVR_enavg       = 0.0
       if (associated(ALBVF_enavg))       ALBVF_enavg       = 0.0
       if (associated(ALBNR_enavg))       ALBNR_enavg       = 0.0
       if (associated(ALBNF_enavg))       ALBNF_enavg       = 0.0
       if (associated(TST_enavg))         TST_enavg         = 0.0
       if (associated(LST_enavg))         LST_enavg         = 0.0
       if (associated(QST_enavg))         QST_enavg         = 0.0
       if (associated(TH_enavg))          TH_enavg          = 0.0
       if (associated(QH_enavg))          QH_enavg          = 0.0
       if (associated(DELTS_enavg))       DELTS_enavg       = 0.0
       if (associated(DELQS_enavg))       DELQS_enavg       = 0.0
       if (associated(CHT_enavg))         CHT_enavg         = 0.0
       if (associated(CMT_enavg))         CMT_enavg         = 0.0
       if (associated(CQT_enavg))         CQT_enavg         = 0.0
       if (associated(CNT_enavg))         CNT_enavg         = 0.0
       if (associated(RIT_enavg))         RIT_enavg         = 0.0
       if (associated(ACCUM_enavg))       ACCUM_enavg       = 0.0
       if (associated(EVPICE_GL_enavg))   EVPICE_GL_enavg   = 0.0
       if (associated(SUBLIM_enavg))      SUBLIM_enavg      = 0.0
       if (associated(SNOMAS_GL_enavg))   SNOMAS_GL_enavg   = 0.0
       if (associated(SNOWMASS_enavg))    SNOWMASS_enavg    = 0.0
       if (associated(SNOWDP_GL_enavg))   SNOWDP_GL_enavg   = 0.0
       if (associated(ASNOW_GL_enavg))    ASNOW_GL_enavg    = 0.0
       if (associated(WESNEXT_enavg))     WESNEXT_enavg     = 0.0
       if (associated(WESNSC_enavg))      WESNSC_enavg      = 0.0
       if (associated(SNDZSC_enavg))      SNDZSC_enavg      = 0.0
       if (associated(WESNPREC_enavg))    WESNPREC_enavg    = 0.0
       if (associated(SNDZPREC_enavg))    SNDZPREC_enavg    = 0.0
       if (associated(SNDZ1PERC_enavg))   SNDZ1PERC_enavg   = 0.0
       if (associated(WESNBOT_enavg))     WESNBOT_enavg     = 0.0
       if (associated(RAINRFZ_enavg))     RAINRFZ_enavg     = 0.0
       if (associated(SMELT_enavg))       SMELT_enavg       = 0.0
       if (associated(IMELT_enavg))       IMELT_enavg       = 0.0
       if (associated(SNOWALB_enavg))     SNOWALB_enavg     = 0.0
       if (associated(SNICEALB_enavg))    SNICEALB_enavg    = 0.0
       if (associated(MELTWTR_enavg))     MELTWTR_enavg     = 0.0
       if (associated(MELTWTRCONT_enavg)) MELTWTRCONT_enavg = 0.0
       if (associated(LWC_enavg))         LWC_enavg         = 0.0
       if (associated(RUNOFF_enavg))      RUNOFF_enavg      = 0.0
       if (associated(GUST_enavg))        GUST_enavg        = 0.0
       if (associated(VENT_enavg))        VENT_enavg        = 0.0
       if (associated(Z0_enavg))          Z0_enavg          = 0.0
       if (associated(Z0H_enavg))         Z0H_enavg         = 0.0
       if (associated(MOT2M_enavg))       MOT2M_enavg       = 0.0
       if (associated(MOQ2M_enavg))       MOQ2M_enavg       = 0.0
       if (associated(MOU2M_enavg))       MOU2M_enavg       = 0.0
       if (associated(MOV2M_enavg))       MOV2M_enavg       = 0.0
       if (associated(MOT10M_enavg))      MOT10M_enavg      = 0.0
       if (associated(MOQ10M_enavg))      MOQ10M_enavg      = 0.0
       if (associated(MOU10M_enavg))      MOU10M_enavg      = 0.0
       if (associated(MOV10M_enavg))      MOV10M_enavg      = 0.0
       if (associated(MOU50M_enavg))      MOU50M_enavg      = 0.0
       if (associated(MOV50M_enavg))      MOV50M_enavg      = 0.0
       if (associated(EVAPOUT_enavg))     EVAPOUT_enavg     = 0.0
       if (associated(SHOUT_enavg))       SHOUT_enavg       = 0.0
       if (associated(HLWUP_enavg))       HLWUP_enavg       = 0.0
       if (associated(LWNDSRF_enavg))     LWNDSRF_enavg     = 0.0
       if (associated(SWNDSRF_enavg))     SWNDSRF_enavg     = 0.0
       if (associated(HLATN_enavg))       HLATN_enavg       = 0.0
       if (associated(DNICFLX_enavg))     DNICFLX_enavg     = 0.0
       if (associated(GHSNOW_enavg))      GHSNOW_enavg      = 0.0
       if (associated(GHTSKIN_enavg))     GHTSKIN_enavg     = 0.0
       if (associated(ITY_enavg))         ITY_enavg         = 0.0
       if (associated(RMELTDU001_enavg))  RMELTDU001_enavg  = 0.0
       if (associated(RMELTDU002_enavg))  RMELTDU002_enavg  = 0.0
       if (associated(RMELTDU003_enavg))  RMELTDU003_enavg  = 0.0
       if (associated(RMELTDU004_enavg))  RMELTDU004_enavg  = 0.0
       if (associated(RMELTDU005_enavg))  RMELTDU005_enavg  = 0.0
       if (associated(RMELTBC001_enavg))  RMELTBC001_enavg  = 0.0
       if (associated(RMELTBC002_enavg))  RMELTBC002_enavg  = 0.0
       if (associated(RMELTOC001_enavg))  RMELTOC001_enavg  = 0.0
       if (associated(RMELTOC002_enavg))  RMELTOC002_enavg  = 0.0
       if (associated(RHOSNOW_enavg))     RHOSNOW_enavg     = 0.0
       if (associated(TSNOW_enavg))       TSNOW_enavg       = 0.0
       if (associated(TICE0_enavg))       TICE0_enavg       = 0.0
       if (associated(WSNOW_enavg))       WSNOW_enavg       = 0.0
       if (associated(ZSNOW_enavg))       ZSNOW_enavg       = 0.0
       if (associated(DRHOS0_enavg))      DRHOS0_enavg      = 0.0
       if (associated(WESNEX_enavg))      WESNEX_enavg      = 0.0
       if (associated(WESNPERC_enavg))    WESNPERC_enavg    = 0.0
       if (associated(WESNDENS_enavg))    WESNDENS_enavg    = 0.0
       if (associated(WESNREPAR_enavg))   WESNREPAR_enavg   = 0.0
    endif

    ! Accumulate ensemble members (1d fields)
    if (associated(EMIS_enavg))        EMIS_enavg        = EMIS_enavg        + EMIS
    if (associated(ALBVR_enavg))       ALBVR_enavg       = ALBVR_enavg       + ALBVR
    if (associated(ALBVF_enavg))       ALBVF_enavg       = ALBVF_enavg       + ALBVF
    if (associated(ALBNR_enavg))       ALBNR_enavg       = ALBNR_enavg       + ALBNR
    if (associated(ALBNF_enavg))       ALBNF_enavg       = ALBNF_enavg       + ALBNF
    if (associated(TST_enavg))         TST_enavg         = TST_enavg         + TST
    if (associated(LST_enavg))         LST_enavg         = LST_enavg         + LST
    if (associated(QST_enavg))         QST_enavg         = QST_enavg         + QST
    if (associated(TH_enavg))          TH_enavg          = TH_enavg          + TH
    if (associated(QH_enavg))          QH_enavg          = QH_enavg          + QH
    if (associated(DELTS_enavg))       DELTS_enavg       = DELTS_enavg       + DELTS
    if (associated(DELQS_enavg))       DELQS_enavg       = DELQS_enavg       + DELQS
    if (associated(CHT_enavg))         CHT_enavg         = CHT_enavg         + CHT
    if (associated(CMT_enavg))         CMT_enavg         = CMT_enavg         + CMT
    if (associated(CQT_enavg))         CQT_enavg         = CQT_enavg         + CQT
    if (associated(CNT_enavg))         CNT_enavg         = CNT_enavg         + CNT
    if (associated(RIT_enavg))         RIT_enavg         = RIT_enavg         + RIT
    if (associated(ACCUM_enavg))       ACCUM_enavg       = ACCUM_enavg       + ACCUM
    if (associated(EVPICE_GL_enavg))   EVPICE_GL_enavg   = EVPICE_GL_enavg   + EVPICE_GL
    if (associated(SUBLIM_enavg))      SUBLIM_enavg      = SUBLIM_enavg      + SUBLIM
    if (associated(SNOMAS_GL_enavg))   SNOMAS_GL_enavg   = SNOMAS_GL_enavg   + SNOMAS_GL
    if (associated(SNOWMASS_enavg))    SNOWMASS_enavg    = SNOWMASS_enavg    + SNOWMASS
    if (associated(SNOWDP_GL_enavg))   SNOWDP_GL_enavg   = SNOWDP_GL_enavg   + SNOWDP_GL
    if (associated(ASNOW_GL_enavg))    ASNOW_GL_enavg    = ASNOW_GL_enavg    + ASNOW_GL
    if (associated(WESNEXT_enavg))     WESNEXT_enavg     = WESNEXT_enavg     + WESNEXT
    if (associated(WESNSC_enavg))      WESNSC_enavg      = WESNSC_enavg      + WESNSC
    if (associated(SNDZSC_enavg))      SNDZSC_enavg      = SNDZSC_enavg      + SNDZSC
    if (associated(WESNPREC_enavg))    WESNPREC_enavg    = WESNPREC_enavg    + WESNPREC
    if (associated(SNDZPREC_enavg))    SNDZPREC_enavg    = SNDZPREC_enavg    + SNDZPREC
    if (associated(SNDZ1PERC_enavg))   SNDZ1PERC_enavg   = SNDZ1PERC_enavg   + SNDZ1PERC
    if (associated(WESNBOT_enavg))     WESNBOT_enavg     = WESNBOT_enavg     + WESNBOT
    if (associated(RAINRFZ_enavg))     RAINRFZ_enavg     = RAINRFZ_enavg     + RAINRFZ
    if (associated(SMELT_enavg))       SMELT_enavg       = SMELT_enavg       + SMELT
    if (associated(IMELT_enavg))       IMELT_enavg       = IMELT_enavg       + IMELT
    if (associated(SNOWALB_enavg))     SNOWALB_enavg     = SNOWALB_enavg     + SNOWALB
    if (associated(SNICEALB_enavg))    SNICEALB_enavg    = SNICEALB_enavg    + SNICEALB
    if (associated(MELTWTR_enavg))     MELTWTR_enavg     = MELTWTR_enavg     + MELTWTR
    if (associated(MELTWTRCONT_enavg)) MELTWTRCONT_enavg = MELTWTRCONT_enavg + MELTWTRCONT
    if (associated(LWC_enavg))         LWC_enavg         = LWC_enavg         + LWC
    if (associated(RUNOFF_enavg))      RUNOFF_enavg      = RUNOFF_enavg      + RUNOFF
    if (associated(GUST_enavg))        GUST_enavg        = GUST_enavg        + GUST
    if (associated(VENT_enavg))        VENT_enavg        = VENT_enavg        + VENT
    if (associated(Z0_enavg))          Z0_enavg          = Z0_enavg          + Z0
    if (associated(Z0H_enavg))         Z0H_enavg         = Z0H_enavg         + Z0H
    if (associated(MOT2M_enavg))       MOT2M_enavg       = MOT2M_enavg       + MOT2M
    if (associated(MOQ2M_enavg))       MOQ2M_enavg       = MOQ2M_enavg       + MOQ2M
    if (associated(MOU2M_enavg))       MOU2M_enavg       = MOU2M_enavg       + MOU2M
    if (associated(MOV2M_enavg))       MOV2M_enavg       = MOV2M_enavg       + MOV2M
    if (associated(MOT10M_enavg))      MOT10M_enavg      = MOT10M_enavg      + MOT10M
    if (associated(MOQ10M_enavg))      MOQ10M_enavg      = MOQ10M_enavg      + MOQ10M
    if (associated(MOU10M_enavg))      MOU10M_enavg      = MOU10M_enavg      + MOU10M
    if (associated(MOV10M_enavg))      MOV10M_enavg      = MOV10M_enavg      + MOV10M
    if (associated(MOU50M_enavg))      MOU50M_enavg      = MOU50M_enavg      + MOU50M
    if (associated(MOV50M_enavg))      MOV50M_enavg      = MOV50M_enavg      + MOV50M
    if (associated(EVAPOUT_enavg))     EVAPOUT_enavg     = EVAPOUT_enavg     + EVAPOUT
    if (associated(SHOUT_enavg))       SHOUT_enavg       = SHOUT_enavg       + SHOUT
    if (associated(HLWUP_enavg))       HLWUP_enavg       = HLWUP_enavg       + HLWUP
    if (associated(LWNDSRF_enavg))     LWNDSRF_enavg     = LWNDSRF_enavg     + LWNDSRF
    if (associated(SWNDSRF_enavg))     SWNDSRF_enavg     = SWNDSRF_enavg     + SWNDSRF
    if (associated(HLATN_enavg))       HLATN_enavg       = HLATN_enavg       + HLATN
    if (associated(DNICFLX_enavg))     DNICFLX_enavg     = DNICFLX_enavg     + DNICFLX
    if (associated(GHSNOW_enavg))      GHSNOW_enavg      = GHSNOW_enavg      + GHSNOW
    if (associated(GHTSKIN_enavg))     GHTSKIN_enavg     = GHTSKIN_enavg     + GHTSKIN
    if (associated(ITY_enavg))         ITY_enavg         = ITY_enavg         + ITY
    if (associated(RMELTDU001_enavg))  RMELTDU001_enavg  = RMELTDU001_enavg  + RMELTDU001
    if (associated(RMELTDU002_enavg))  RMELTDU002_enavg  = RMELTDU002_enavg  + RMELTDU002
    if (associated(RMELTDU003_enavg))  RMELTDU003_enavg  = RMELTDU003_enavg  + RMELTDU003
    if (associated(RMELTDU004_enavg))  RMELTDU004_enavg  = RMELTDU004_enavg  + RMELTDU004
    if (associated(RMELTDU005_enavg))  RMELTDU005_enavg  = RMELTDU005_enavg  + RMELTDU005
    if (associated(RMELTBC001_enavg))  RMELTBC001_enavg  = RMELTBC001_enavg  + RMELTBC001
    if (associated(RMELTBC002_enavg))  RMELTBC002_enavg  = RMELTBC002_enavg  + RMELTBC002
    if (associated(RMELTOC001_enavg))  RMELTOC001_enavg  = RMELTOC001_enavg  + RMELTOC001
    if (associated(RMELTOC002_enavg))  RMELTOC002_enavg  = RMELTOC002_enavg  + RMELTOC002

    ! Accumulate ensemble members (multi-dimensional fields)
    if (associated(RHOSNOW_enavg))     RHOSNOW_enavg     = RHOSNOW_enavg     + RHOSNOW
    if (associated(TSNOW_enavg))       TSNOW_enavg       = TSNOW_enavg       + TSNOW
    if (associated(TICE0_enavg))       TICE0_enavg       = TICE0_enavg       + TICE0
    if (associated(WSNOW_enavg))       WSNOW_enavg       = WSNOW_enavg       + WSNOW
    if (associated(ZSNOW_enavg))       ZSNOW_enavg       = ZSNOW_enavg       + ZSNOW
    if (associated(DRHOS0_enavg))      DRHOS0_enavg      = DRHOS0_enavg      + DRHOS0
    if (associated(WESNEX_enavg))      WESNEX_enavg      = WESNEX_enavg      + WESNEX
    if (associated(WESNPERC_enavg))    WESNPERC_enavg    = WESNPERC_enavg    + WESNPERC
    if (associated(WESNDENS_enavg))    WESNDENS_enavg    = WESNDENS_enavg    + WESNDENS
    if (associated(WESNREPAR_enavg))   WESNREPAR_enavg   = WESNREPAR_enavg   + WESNREPAR

    collect_landice_counter = collect_landice_counter + 1

    if (collect_landice_counter == NUM_ENSEMBLE) then
       collect_landice_counter = 0

       ! Divide by NUM_ENSEMBLE (1d fields)
       if (associated(EMIS_enavg))        EMIS_enavg        = EMIS_enavg        / NUM_ENSEMBLE
       if (associated(ALBVR_enavg))       ALBVR_enavg       = ALBVR_enavg       / NUM_ENSEMBLE
       if (associated(ALBVF_enavg))       ALBVF_enavg       = ALBVF_enavg       / NUM_ENSEMBLE
       if (associated(ALBNR_enavg))       ALBNR_enavg       = ALBNR_enavg       / NUM_ENSEMBLE
       if (associated(ALBNF_enavg))       ALBNF_enavg       = ALBNF_enavg       / NUM_ENSEMBLE
       if (associated(TST_enavg))         TST_enavg         = TST_enavg         / NUM_ENSEMBLE
       if (associated(LST_enavg))         LST_enavg         = LST_enavg         / NUM_ENSEMBLE
       if (associated(QST_enavg))         QST_enavg         = QST_enavg         / NUM_ENSEMBLE
       if (associated(TH_enavg))          TH_enavg          = TH_enavg          / NUM_ENSEMBLE
       if (associated(QH_enavg))          QH_enavg          = QH_enavg          / NUM_ENSEMBLE
       if (associated(DELTS_enavg))       DELTS_enavg       = DELTS_enavg       / NUM_ENSEMBLE
       if (associated(DELQS_enavg))       DELQS_enavg       = DELQS_enavg       / NUM_ENSEMBLE
       if (associated(CHT_enavg))         CHT_enavg         = CHT_enavg         / NUM_ENSEMBLE
       if (associated(CMT_enavg))         CMT_enavg         = CMT_enavg         / NUM_ENSEMBLE
       if (associated(CQT_enavg))         CQT_enavg         = CQT_enavg         / NUM_ENSEMBLE
       if (associated(CNT_enavg))         CNT_enavg         = CNT_enavg         / NUM_ENSEMBLE
       if (associated(RIT_enavg))         RIT_enavg         = RIT_enavg         / NUM_ENSEMBLE
       if (associated(ACCUM_enavg))       ACCUM_enavg       = ACCUM_enavg       / NUM_ENSEMBLE
       if (associated(EVPICE_GL_enavg))   EVPICE_GL_enavg   = EVPICE_GL_enavg   / NUM_ENSEMBLE
       if (associated(SUBLIM_enavg))      SUBLIM_enavg      = SUBLIM_enavg      / NUM_ENSEMBLE
       if (associated(SNOMAS_GL_enavg))   SNOMAS_GL_enavg   = SNOMAS_GL_enavg   / NUM_ENSEMBLE
       if (associated(SNOWMASS_enavg))    SNOWMASS_enavg    = SNOWMASS_enavg    / NUM_ENSEMBLE
       if (associated(SNOWDP_GL_enavg))   SNOWDP_GL_enavg   = SNOWDP_GL_enavg   / NUM_ENSEMBLE
       if (associated(ASNOW_GL_enavg))    ASNOW_GL_enavg    = ASNOW_GL_enavg    / NUM_ENSEMBLE
       if (associated(WESNEXT_enavg))     WESNEXT_enavg     = WESNEXT_enavg     / NUM_ENSEMBLE
       if (associated(WESNSC_enavg))      WESNSC_enavg      = WESNSC_enavg      / NUM_ENSEMBLE
       if (associated(SNDZSC_enavg))      SNDZSC_enavg      = SNDZSC_enavg      / NUM_ENSEMBLE
       if (associated(WESNPREC_enavg))    WESNPREC_enavg    = WESNPREC_enavg    / NUM_ENSEMBLE
       if (associated(SNDZPREC_enavg))    SNDZPREC_enavg    = SNDZPREC_enavg    / NUM_ENSEMBLE
       if (associated(SNDZ1PERC_enavg))   SNDZ1PERC_enavg   = SNDZ1PERC_enavg   / NUM_ENSEMBLE
       if (associated(WESNBOT_enavg))     WESNBOT_enavg     = WESNBOT_enavg     / NUM_ENSEMBLE
       if (associated(RAINRFZ_enavg))     RAINRFZ_enavg     = RAINRFZ_enavg     / NUM_ENSEMBLE
       if (associated(SMELT_enavg))       SMELT_enavg       = SMELT_enavg       / NUM_ENSEMBLE
       if (associated(IMELT_enavg))       IMELT_enavg       = IMELT_enavg       / NUM_ENSEMBLE
       if (associated(SNOWALB_enavg))     SNOWALB_enavg     = SNOWALB_enavg     / NUM_ENSEMBLE
       if (associated(SNICEALB_enavg))    SNICEALB_enavg    = SNICEALB_enavg    / NUM_ENSEMBLE
       if (associated(MELTWTR_enavg))     MELTWTR_enavg     = MELTWTR_enavg     / NUM_ENSEMBLE
       if (associated(MELTWTRCONT_enavg)) MELTWTRCONT_enavg = MELTWTRCONT_enavg / NUM_ENSEMBLE
       if (associated(LWC_enavg))         LWC_enavg         = LWC_enavg         / NUM_ENSEMBLE
       if (associated(RUNOFF_enavg))      RUNOFF_enavg      = RUNOFF_enavg      / NUM_ENSEMBLE
       if (associated(GUST_enavg))        GUST_enavg        = GUST_enavg        / NUM_ENSEMBLE
       if (associated(VENT_enavg))        VENT_enavg        = VENT_enavg        / NUM_ENSEMBLE
       if (associated(Z0_enavg))          Z0_enavg          = Z0_enavg          / NUM_ENSEMBLE
       if (associated(Z0H_enavg))         Z0H_enavg         = Z0H_enavg         / NUM_ENSEMBLE
       if (associated(MOT2M_enavg))       MOT2M_enavg       = MOT2M_enavg       / NUM_ENSEMBLE
       if (associated(MOQ2M_enavg))       MOQ2M_enavg       = MOQ2M_enavg       / NUM_ENSEMBLE
       if (associated(MOU2M_enavg))       MOU2M_enavg       = MOU2M_enavg       / NUM_ENSEMBLE
       if (associated(MOV2M_enavg))       MOV2M_enavg       = MOV2M_enavg       / NUM_ENSEMBLE
       if (associated(MOT10M_enavg))      MOT10M_enavg      = MOT10M_enavg      / NUM_ENSEMBLE
       if (associated(MOQ10M_enavg))      MOQ10M_enavg      = MOQ10M_enavg      / NUM_ENSEMBLE
       if (associated(MOU10M_enavg))      MOU10M_enavg      = MOU10M_enavg      / NUM_ENSEMBLE
       if (associated(MOV10M_enavg))      MOV10M_enavg      = MOV10M_enavg      / NUM_ENSEMBLE
       if (associated(MOU50M_enavg))      MOU50M_enavg      = MOU50M_enavg      / NUM_ENSEMBLE
       if (associated(MOV50M_enavg))      MOV50M_enavg      = MOV50M_enavg      / NUM_ENSEMBLE
       if (associated(EVAPOUT_enavg))     EVAPOUT_enavg     = EVAPOUT_enavg     / NUM_ENSEMBLE
       if (associated(SHOUT_enavg))       SHOUT_enavg       = SHOUT_enavg       / NUM_ENSEMBLE
       if (associated(HLWUP_enavg))       HLWUP_enavg       = HLWUP_enavg       / NUM_ENSEMBLE
       if (associated(LWNDSRF_enavg))     LWNDSRF_enavg     = LWNDSRF_enavg     / NUM_ENSEMBLE
       if (associated(SWNDSRF_enavg))     SWNDSRF_enavg     = SWNDSRF_enavg     / NUM_ENSEMBLE
       if (associated(HLATN_enavg))       HLATN_enavg       = HLATN_enavg       / NUM_ENSEMBLE
       if (associated(DNICFLX_enavg))     DNICFLX_enavg     = DNICFLX_enavg     / NUM_ENSEMBLE
       if (associated(GHSNOW_enavg))      GHSNOW_enavg      = GHSNOW_enavg      / NUM_ENSEMBLE
       if (associated(GHTSKIN_enavg))     GHTSKIN_enavg     = GHTSKIN_enavg     / NUM_ENSEMBLE
       if (associated(ITY_enavg))         ITY_enavg         = ITY_enavg         / NUM_ENSEMBLE
       if (associated(RMELTDU001_enavg))  RMELTDU001_enavg  = RMELTDU001_enavg  / NUM_ENSEMBLE
       if (associated(RMELTDU002_enavg))  RMELTDU002_enavg  = RMELTDU002_enavg  / NUM_ENSEMBLE
       if (associated(RMELTDU003_enavg))  RMELTDU003_enavg  = RMELTDU003_enavg  / NUM_ENSEMBLE
       if (associated(RMELTDU004_enavg))  RMELTDU004_enavg  = RMELTDU004_enavg  / NUM_ENSEMBLE
       if (associated(RMELTDU005_enavg))  RMELTDU005_enavg  = RMELTDU005_enavg  / NUM_ENSEMBLE
       if (associated(RMELTBC001_enavg))  RMELTBC001_enavg  = RMELTBC001_enavg  / NUM_ENSEMBLE
       if (associated(RMELTBC002_enavg))  RMELTBC002_enavg  = RMELTBC002_enavg  / NUM_ENSEMBLE
       if (associated(RMELTOC001_enavg))  RMELTOC001_enavg  = RMELTOC001_enavg  / NUM_ENSEMBLE
       if (associated(RMELTOC002_enavg))  RMELTOC002_enavg  = RMELTOC002_enavg  / NUM_ENSEMBLE

       ! Divide by NUM_ENSEMBLE (multi-dimensional fields)
       if (associated(RHOSNOW_enavg))     RHOSNOW_enavg     = RHOSNOW_enavg     / NUM_ENSEMBLE
       if (associated(TSNOW_enavg))       TSNOW_enavg       = TSNOW_enavg       / NUM_ENSEMBLE
       if (associated(TICE0_enavg))       TICE0_enavg       = TICE0_enavg       / NUM_ENSEMBLE
       if (associated(WSNOW_enavg))       WSNOW_enavg       = WSNOW_enavg       / NUM_ENSEMBLE
       if (associated(ZSNOW_enavg))       ZSNOW_enavg       = ZSNOW_enavg       / NUM_ENSEMBLE
       if (associated(DRHOS0_enavg))      DRHOS0_enavg      = DRHOS0_enavg      / NUM_ENSEMBLE
       if (associated(WESNEX_enavg))      WESNEX_enavg      = WESNEX_enavg      / NUM_ENSEMBLE
       if (associated(WESNPERC_enavg))    WESNPERC_enavg    = WESNPERC_enavg    / NUM_ENSEMBLE
       if (associated(WESNDENS_enavg))    WESNDENS_enavg    = WESNDENS_enavg    / NUM_ENSEMBLE
       if (associated(WESNREPAR_enavg))   WESNREPAR_enavg   = WESNREPAR_enavg   / NUM_ENSEMBLE
    endif

    call MAPL_TimerOff(MAPL, "Collect_landice")
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

end module GEOS_LandiceAvgGridCompMod
