esma_set_this ()

set (alldirs
  GEOSmetforce_GridComp
  GEOSlandpert_GridComp
  GEOSens_GridComp
  GEOSlandassim_GridComp
  )

esma_add_library(${this}
  SRCS GEOS_LdasGridComp.F90
  SUBCOMPONENTS ${alldirs}
  SUBDIRS LDAS_Shared
  DEPENDENCIES GEOSland_GridComp GEOSlandice_GridComp makebcs MAPL
  INCLUDES ${INC_ESMF})

esma_add_subdirectory(GEOSldas_App)
