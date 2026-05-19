esma_set_this ()

set (SRCS
   GEOS_EnsGridComp.F90
   GEOS_ForceAvgGridComp.F90
   GEOS_LandiceAvgGridComp.F90
   GEOS_LandAvgGridComp.F90
   GEOS_RouteAvgGridComp.F90
)

esma_add_library(${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSland_GridComp GEOSlandice_GridComp GEOSroute_GridComp GEOS_LdasShared MAPL
  INCLUDES ${INC_ESMF})
