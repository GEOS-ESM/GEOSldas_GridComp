# -----------------------------------------------------------------------
# compare_pert/CMakeLists.txt
#
# Builds the MPI comparison program for flat lat-lon FFT and sphere-native
# SMERFS perturbations on a 360x180 grid.
#
# Executable:
#   compare_pert.x      - runs both perturbation paths and writes binary output
#
# Python script (installed to ${CMAKE_INSTALL_PREFIX}/bin):
#   compare_pert.py     - reads both .bin outputs, produces comparison figure
#
# Both executables link against GEOSlandpert_GridComp, which transitively
# carries MAPL, ESMF, MKL, smerfs_f, MPI, netCDF, and HDF5.
# The MKL_AVAILABLE compile definition must be repeated here because the
# parent library applies it only as PRIVATE.
# -----------------------------------------------------------------------

# ---- Combined comparison executable ------------------------------------
ecbuild_add_executable(
  TARGET  compare_pert.x
  SOURCES compare_pert.F90
  LIBS    GEOSlandpert_GridComp MAPL)

ecbuild_add_executable(
  TARGET  pert_statistics.x
  SOURCES pert_statistics.F90
  LIBS    GEOSlandpert_GridComp MAPL)

#ecbuild_add_executable(
#  TARGET  compare_timing.x
#  SOURCES compare_timing.F90
#  LIBS    GEOSlandpert_GridComp MAPL)

# ---- MKL compile definition (matches parent library treatment) ---------
if (MKL_FOUND AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_definitions(compare_pert.x PRIVATE MKL_AVAILABLE)
  target_compile_definitions(pert_statistics.x PRIVATE MKL_AVAILABLE)
#  target_compile_definitions(compare_timing.x PRIVATE MKL_AVAILABLE)
endif ()

# ---- Install Python comparison script ----------------------------------
install(
  PROGRAMS compare_pert.py
  DESTINATION bin)
