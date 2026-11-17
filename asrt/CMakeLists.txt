cmake_minimum_required(VERSION 3.12)
project(asrt LANGUAGES Fortran)

# ASRT: Active Sensor Radiative Transfer, a radar reflectivity operator (no CRTM dependency).
# Needs netCDF-Fortran only to read bulk scattering lookup tables.

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(PkgConfig QUIET)
find_path(NETCDF_F90_INCLUDE netcdf.mod HINTS $ENV{NETCDF}/include $ENV{NETCDF_FORTRAN}/include
          /usr/include /usr/local/include)
find_library(NETCDF_F90_LIB netcdff HINTS $ENV{NETCDF}/lib $ENV{NETCDF_FORTRAN}/lib)
find_library(NETCDF_C_LIB netcdf HINTS $ENV{NETCDF}/lib $ENV{NETCDF_C}/lib)

set(ASRT_SRC
  src/asrt_kinds.f90
  src/asrt_interp.f90
  src/asrt_gas.f90
  src/asrt_hydro_approx.f90
  src/asrt_species.f90
  src/asrt_tdts.F90
  src/asrt_operator.f90
  src/asrt_simple.f90)

add_library(asrt ${ASRT_SRC})
set_target_properties(asrt PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module)
target_include_directories(asrt PUBLIC ${CMAKE_BINARY_DIR}/module ${NETCDF_F90_INCLUDE})
target_link_libraries(asrt PUBLIC ${NETCDF_F90_LIB} ${NETCDF_C_LIB})

option(ASRT_BUILD_TESTS "Build the tests" ON)
if(ASRT_BUILD_TESTS)
  enable_testing()
  add_executable(test_asrt test/test_asrt.f90)
  target_link_libraries(test_asrt asrt)
  add_test(NAME asrt_tl_ad_k COMMAND test_asrt ${ASRT_CLOUDCOEFF})
endif()
