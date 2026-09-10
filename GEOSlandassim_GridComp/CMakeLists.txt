esma_set_this ()

set (alldirs
  GEOSexportcatchincr_GridComp 
  )

set (SRCS
  hdf4_c_api.c hdf4_fortran_api.F90 io_hdf5.F90 enkf_general.F90 adapt_types.F90 mwRTM_types.F90 clsm_ensupd_glob_param.F90
  mwRTM_routines.F90 clsm_ensupd_upd_routines.F90 clsm_ensdrv_drv_routines.F90
  clsm_ensupd_read_obs.F90 catch_bias_types.F90 clsm_bias_routines.F90 clsm_adapt_routines.F90
  clsm_ensupd_enkf_update.F90 clsm_ensdrv_out_routines.F90 GEOS_LandAssimGridComp.F90
  )

find_package(HDF5 REQUIRED COMPONENTS Fortran)
find_path(HDF4_INCLUDE_DIR hdf.h PATH_SUFFIXES hdf REQUIRED)
find_library(HDF4_MFHDF_LIBRARY NAMES mfhdf REQUIRED)
find_library(HDF4_DF_LIBRARY NAMES df REQUIRED)

esma_add_library (${this}
  SRCS ${SRCS}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES GEOS_LdasShared GEOSensavg_GridComp GEOSlandpert_GridComp GEOSland_GridComp makebcs MAPL NCEP_bufr_r4i4 GMAO_gfio_r4 hdf5hl_fortran hdf5_fortran ${HDF4_MFHDF_LIBRARY} ${HDF4_DF_LIBRARY} ${NETCDF_LIBRARIES}
  INCLUDES ${INC_ESMF} ${INC_HDF5} ${HDF4_INCLUDE_DIR})

target_compile_definitions (${this} PRIVATE LDAS_MPI)

# NAG: HDF4 Fortran API (VSFREAD) is called with varying types intentionally
if (CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_options(${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-wmismatch=VSFREAD>)
endif()
