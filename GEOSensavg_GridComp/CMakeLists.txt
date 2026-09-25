esma_set_this ()

set (SRCS
   GEOS_MetforceAvgGridComp.F90
   GEOS_LandiceAvgGridComp.F90
   GEOS_LakeAvgGridComp.F90
   GEOS_LandAvgGridComp.F90		
   GEOS_RouteAvgGridComp.F90
)

esma_add_library(${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSland_GridComp GEOSlake_GridComp GEOSlandice_GridComp GEOSroute_GridComp GEOS_LdasShared MAPL
  INCLUDES ${INC_ESMF})
