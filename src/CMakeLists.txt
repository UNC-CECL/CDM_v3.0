cmake_minimum_required(VERSION 3.0)

project (cdm CXX)

set (CDM_VERSION 3.0)

find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
include_directories(${FFTW3_INCLUDE_DIRS})
link_directories(${FFTW3_LIBRARY_DIRS})

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
target_link_libraries(coastal-dune-model fftw3)

install(
  TARGETS coastal-dune-model
  RUNTIME DESTINATION bin
)
