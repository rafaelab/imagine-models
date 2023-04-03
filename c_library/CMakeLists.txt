cmake_minimum_required(VERSION 3.4...3.18)

project(ImagineModels VERSION 0.0.1 DESCRIPTION "A library for Galactic models")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/findFFTW")
include(GNUInstallDirs)

FIND_PACKAGE(FFTW)
FIND_PACKAGE(autodiff)


set(INCLUDE_INSTALL_DIR
"${CMAKE_INSTALL_INCLUDEDIR}/ImagineModels"
# CACHE PATH "The directory relative to CMAKE_INSTALL_PREFIX where ImagineModel header files are installed"
)

if (FFTW_FOUND)
    set(IM_PUBLIC_HEADERS
        "headers/RegularModels.h"
        "headers/RandomModels.h"

        "headers/exceptions.h"

        "headers/Field.h"
        "headers/RegularField.h"
        "headers/RandomField.h"
        "headers/random.tpp"

        "headers/helpers.h"
        "headers/Uniform.h"
        "headers/RegularJF12.h"
        "headers/Jaffe.h"
        "headers/Helix.h"
        "headers/YMW.h"
        "headers/Sun.h"
        "headers/TinyakovTkachev.h"
        "headers/HarariMollerachRoulet.h"

        "headers/GaussianScalar.h"
        "headers/LogNormal.h"
        "headers/RandomJF12.h"
        "headers/EnsslinSteininger.h"
    )
    add_library(ImagineModels SHARED
        source/helix.cc
        source/jaffe.cc
        source/regularjf12.cc
        source/gaussianscalar.cc
        source/lognormal.cc
        source/randomjf12.cc
        source/ensslinsteininger.cc
        source/ymw16.cc
        source/tinyakovtkachev.cc
        source/hararimollerachroulet.cc
        source/sun.cc

    )
else()
    set(IM_PUBLIC_HEADERS
        "headers/RegularModels.h"
        "headers/exceptions.h"
        "headers/Field.h"
        "headers/RegularField.h"

        "headers/helpers.h"
        "headers/Uniform.h"
        "headers/RegularJF12.h"
        "headers/Jaffe.h"
        "headers/Helix.h"
        "headers/YMW.h"
        "headers/TinyakovTkachev.h"
        "headers/Sun.h"
        "headers/StanevBSS.h"
        "headers/HarariMollerachRoulet.h"
    )
    add_library(ImagineModels SHARED
    source/helix.cc
    source/jaffe.cc
    source/regularjf12.cc
    source/ymw16.cc
    source/sun.cc
    source/stanevbss.cc
    source/tinyakovtkachev.cc
    source/hararimollerachroulet.cc
    )
endif()

target_include_directories(ImagineModels PRIVATE .)

message(FFTW_FOUND="${FFTW_FOUND}")
message(PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
message(IM_PUBLIC_HEADERS="${IM_PUBLIC_HEADERS}")
message(CMAKE_INSTALL_LIBDIR="${CMAKE_INSTALL_LIBDIR}")
message(CMAKE_INSTALL_INCLUDEDIR="${CMAKE_INSTALL_INCLUDEDIR}")
message(CMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
message(CMAKE_MODULE_PATH="${CMAKE_MODULE_PATH}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set_target_properties(ImagineModels PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
    PUBLIC_HEADER "${IM_PUBLIC_HEADERS}")
configure_file(ImagineModels.pc.in ImagineModels.pc @ONLY)
target_compile_features(ImagineModels PRIVATE cxx_std_17)

install(TARGETS ImagineModels
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ImagineModels)
install(FILES ${CMAKE_BINARY_DIR}/ImagineModels.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
