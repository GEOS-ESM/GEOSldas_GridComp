# compare_pert.x: performance and output comparison of rectangle_random_fields
# and sphere_random_fields.  See run_compare.sh.

ecbuild_add_executable(
  TARGET  compare_pert.x
  SOURCES compare_pert.F90
  LIBS    GEOSlandpert_GridComp MAPL)

# MKL_AVAILABLE is PRIVATE in the parent library, so repeat it here.
if (MKL_FOUND AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_definitions(compare_pert.x PRIVATE MKL_AVAILABLE)
endif ()

install(
  PROGRAMS compare_pert.py plot_timing.py run_compare.sh
  DESTINATION bin)
