cmake_minimum_required(VERSION 3.12)
project(radar_operator LANGUAGES Fortran)

# Standalone radar reflectivity operator (no CRTM dependency).
# Needs netCDF-Fortran only to read bulk scattering lookup tables.

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(PkgConfig QUIET)
find_path(NETCDF_F90_INCLUDE netcdf.mod HINTS $ENV{NETCDF}/include $ENV{NETCDF_FORTRAN}/include
          /usr/include /usr/local/include)
find_library(NETCDF_F90_LIB netcdff HINTS $ENV{NETCDF}/lib $ENV{NETCDF_FORTRAN}/lib)
find_library(NETCDF_C_LIB netcdf HINTS $ENV{NETCDF}/lib $ENV{NETCDF_C}/lib)

set(RADAR_SRC
  src/radar_kinds.f90
  src/radar_interp.f90
  src/radar_gas.f90
  src/radar_hydro_approx.f90
  src/radar_species.f90
  src/radar_tdts.F90
  src/radar_operator.f90
  src/radar_simple.f90)

add_library(radar_operator ${RADAR_SRC})
set_target_properties(radar_operator PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module)
target_include_directories(radar_operator PUBLIC ${CMAKE_BINARY_DIR}/module ${NETCDF_F90_INCLUDE})
target_link_libraries(radar_operator PUBLIC ${NETCDF_F90_LIB} ${NETCDF_C_LIB})

option(RADAR_BUILD_TESTS "Build the tests" ON)
if(RADAR_BUILD_TESTS)
  enable_testing()
  add_executable(test_radar_operator test/test_radar_operator.f90)
  target_link_libraries(test_radar_operator radar_operator)
  add_test(NAME radar_operator_tl_ad_k COMMAND test_radar_operator ${RADAR_CLOUDCOEFF})
endif()
