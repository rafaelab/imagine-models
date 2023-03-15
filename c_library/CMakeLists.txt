cmake_minimum_required(VERSION 3.4...3.18)

project(ImagineModels VERSION 0.0.1 DESCRIPTION "A library for Galactic models")

include(GNUInstallDirs)
set(INCLUDE_INSTALL_DIR
"${CMAKE_INSTALL_INCLUDEDIR}/ImagineModels"
# CACHE PATH "The directory relative to CMAKE_INSTALL_PREFIX where ImagineModel header files are installed"
)

set(IM_PUBLIC_HEADERS
"headers/RegularModels.h"
"headers/RandomModels.h"

"headers/exceptions.h"

"headers/Field.h"
"headers/RegularField.h"
"headers/RandomField.h"

"headers/RegularJF12.h"
"headers/Jaffe.h"
"headers/Helix.h"

"headers/RandomJF12.h"
"headers/EnsslinSteininger.h"

    )

add_library(ImagineModels SHARED
#    source/component.cc
#    source/disc.cc
    source/helix.cc
    source/jaffe.cc
    source/regularjf12.cc
    source/randomscalar.cc
    source/randomvector.cc
    source/randomjf12.cc
    source/ensslinsteininger.cc
#    source/ymw16.cc
#    source/regular.cc
)
target_include_directories(ImagineModels PRIVATE .)

message(PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
message(IM_PUBLIC_HEADERS="${IM_PUBLIC_HEADERS}")
message(CMAKE_INSTALL_LIBDIR="${CMAKE_INSTALL_LIBDIR}")
message(CMAKE_INSTALL_INCLUDEDIR="${CMAKE_INSTALL_INCLUDEDIR}")
message(CMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")

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
