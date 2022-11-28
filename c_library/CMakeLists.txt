cmake_minimum_required(VERSION 3.4...3.18)

project(ImagineModels VERSION 0.0.1 DESCRIPTION "A library for parametric Galactic  models")

include(GNUInstallDirs)

add_library(ImagineModels SHARED
#    source/component.cc
#    source/disc.cc
#    source/helix.cc
#    source/jaffe.cc
    source/jf12.cc
#    source/ymw16.cc
#    source/regular.cc
)
set_target_properties(ImagineModels PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
  PUBLIC_HEADER "headers/MagneticField.h;")
configure_file(ImagineModels.pc.in ImagineModels.pc @ONLY)
target_compile_features(ImagineModels PRIVATE cxx_std_17)
target_include_directories(ImagineModels PRIVATE .)
install(TARGETS ImagineModels
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/ImagineModels.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
