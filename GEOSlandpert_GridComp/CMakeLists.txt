esma_set_this ()

set (SRCS
  nr_ran2_gasdev.F90 nr_jacobi.F90 nr_fft.F90
  abstract_random_fields.F90 rectangle_random_fields.F90 sphere_riccati_init.F90 sphere_random_fields.F90
  land_pert.F90 force_and_cat_progn_pert_types.F90 LDAS_PertRoutines.F90
  GEOS_LandPertGridComp.F90
  )

# The Riccati initializer and restartable NR production RNG do not need SMERFS.
# Keep the old MT/Ziggurat seed-only interface available as an opt-in reference.
option (GEOS_SPHERE_ENABLE_SMERFS_REFERENCE
  "Enable the legacy spherical MT/Ziggurat reference RNG (requires @smerfs)" OFF)

set (_sphere_reference_libraries "")
if (GEOS_SPHERE_ENABLE_SMERFS_REFERENCE)
  if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/@smerfs/CMakeLists.txt")
    message (FATAL_ERROR
      "GEOS_SPHERE_ENABLE_SMERFS_REFERENCE=ON requires @smerfs/CMakeLists.txt. "
      "Provide the author's SMERFS checkout or set this option OFF for NR production runs.")
  endif ()
  esma_add_subdirectory(@smerfs)
  if (NOT TARGET smerfs_f)
    message (FATAL_ERROR "The supplied @smerfs checkout did not define the smerfs_f target.")
  endif ()
  list (APPEND _sphere_reference_libraries smerfs_f)
endif ()

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOS_LdasShared GEOSensavg_GridComp GEOSland_GridComp MAPL ${_sphere_reference_libraries} ${MKL_LIBRARIES}
  INCLUDES ${INC_ESMF} ${MKL_INCLUDE_DIRS})

if (GEOS_SPHERE_ENABLE_SMERFS_REFERENCE)
  target_compile_definitions(${this} PRIVATE SPHERE_SMERFS_REFERENCE)
endif ()
unset (_sphere_reference_libraries)

# Enable MKL paths only when MKL was detected and the compiler supports it.
if (MKL_FOUND AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_definitions(${this} PRIVATE MKL_AVAILABLE)
endif ()

esma_add_subdirectory(compare_pert)
