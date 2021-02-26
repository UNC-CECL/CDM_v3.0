cmake_minimum_required(VERSION 3.0)
include(CheckIncludeFile)

# project (cdm CXX)
project (cdm)

set (CDM_VERSION 3.0)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})

check_include_file(unistd.h WITH_UNISTD)
check_include_file(sys/utsname.h WITH_UTSNAME)
find_library(LIB_M m)

if (WITH_UNISTD)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_UNISTD")
endif (WITH_UNISTD)
if (WITH_UTSNAME)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_UTSNAME")
endif (WITH_UTSNAME)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_USE_MATH_DEFINES")

set(
  CDM_SRCS
    PTG_Func2dScalar.cc
    PTG_Func2dVec.cc
    analyze_new.cc
    avalanche.cc
    dune_evolution.cc
    evolution.cc
    flux_stationary.cc
    func.cc
    globals.cc
    influx.cc
    initsurf.cc
    initsurfalea.cc
    initsurfbeach.cc
    iter_aval_new.cc
    rfftw12d.cc
    rotatematrix.cc
    save.cc
    sepbubble.cc
    shear.cc
    shear_hlr.cc
    shore.cc
    storm.cc
    vegetation.cc
    wind.cc
)

add_executable(coastal-dune-model main.cc ${CDM_SRCS})
if (LIB_M)
  target_link_libraries(coastal-dune-model fftw3 m)
elseif (NOT LIB_M)
  target_link_libraries(coastal-dune-model fftw3)
endif (LIB_M)

install(
  TARGETS coastal-dune-model
  RUNTIME DESTINATION bin
)

