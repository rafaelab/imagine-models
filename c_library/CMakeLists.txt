cmake_minimum_required(VERSION 3.4...3.18)

project(ImagineModels VERSION 0.0.1 DESCRIPTION "A library for Galactic models")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/findFFTW")
include(GNUInstallDirs)
include(CTest)

FIND_PACKAGE(FFTW)
FIND_PACKAGE(autodiff)


set(INCLUDE_INSTALL_DIR
"${CMAKE_INSTALL_INCLUDEDIR}/ImagineModels"
# CACHE PATH "The directory relative to CMAKE_INSTALL_PREFIX where ImagineModel header files are installed"
)


set(IM_PUBLIC_HEADERS
    # Generic headers
    "headers/Field.h"
    "headers/RegularField.h"
    "headers/RegularModels.h" # allows the include of all regular models
   
    "headers/exceptions.h"
    "headers/helpers.h"
    "headers/hamunits.h"

    "headers/Uniform.h"
    "headers/RegularJF12.h"
    "headers/Jaffe.h"
    "headers/Helix.h"
    "headers/YMW.h"
    "headers/TinyakovTkachev.h"
    "headers/Sun.h"
    "headers/StanevBSS.h"
    "headers/Fauvet.h"
    "headers/WMAP.h"
    "headers/HarariMollerachRoulet.h"
)

add_library(ImagineModels SHARED
source/helix.cc
source/jaffe.cc
source/regularjf12.cc
source/ymw16.cc
source/sun.cc
source/wmap.cc
source/stanevbss.cc
source/fauvet.cc
source/tinyakovtkachev.cc
source/hararimollerachroulet.cc
)

if (FFTW_FOUND)
    set(_RANDOM_HEADERS
        "headers/RandomModels.h"

        "headers/RandomField.h"
        "headers/random.tpp"

        "headers/GaussianScalar.h"
        "headers/LogNormal.h"
        "headers/RandomJF12.h"
        "headers/EnsslinSteininger.h"
    )
    list(APPEND IM_PUBLIC_HEADERS "${_RANDOM_HEADERS}")
    
    add_library(ImagineModelsRandom SHARED
        source/gaussianscalar.cc
        source/lognormal.cc
        source/randomjf12.cc
        source/ensslinsteininger.cc
    )
    target_link_libraries(ImagineModels PRIVATE ImagineModelsRandom)


endif()

target_include_directories(ImagineModels PRIVATE .)

set(TESTSOURCES grid parameter_update positions)
foreach(test ${TESTSOURCES})
    set(test_name ${test}_ctest)
    add_executable(${test_name} "${PROJECT_SOURCE_DIR}/test/regular/test_${test}.cc")   
    target_link_libraries(${test_name} ImagineModels)
    add_test(${test_name} ${test_name})
endforeach()

enable_testing()


message(TEST_LIBRARY="${PROJECT_SOURCE_DIR}/test/test_regular_fields.cc")
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
