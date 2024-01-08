project(tImagineModels VERSION 0.2 DESCRIPTION "A library for Galactic models.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

set(IM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/headers/)
INCLUDE_DIRECTORIES ( "/usr/local/include/eigen3 " )
set(IM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/source)


set(INCLUDE_INSTALL_DIR
"${CMAKE_INSTALL_INCLUDEDIR}/ImagineModels")
# CACHE PATH "The directory relative to CMAKE_INSTALL_PREFIX where ImagineModel header files are installed"


set(IM_SRC_FILES
  #  ${IM_SOURCE_DIR}/helpers.cc
    ${IM_SOURCE_DIR}/helix.cc
    ${IM_SOURCE_DIR}/jaffe.cc
    ${IM_SOURCE_DIR}/regularjf12.cc
    ${IM_SOURCE_DIR}/ymw16.cc
    ${IM_SOURCE_DIR}/tf17.cc
    ${IM_SOURCE_DIR}/sun.cc
    ${IM_SOURCE_DIR}/pshirkov.cc
    ${IM_SOURCE_DIR}/archimedes.cc
    ${IM_SOURCE_DIR}/han.cc
    ${IM_SOURCE_DIR}/wmap.cc
    ${IM_SOURCE_DIR}/stanevbss.cc
    ${IM_SOURCE_DIR}/fauvet.cc
    ${IM_SOURCE_DIR}/tinyakovtkachev.cc
    ${IM_SOURCE_DIR}/hararimollerachroulet.cc   
)

if(FFTW_FOUND)
    set(IM_SRC_FILES
        ${IM_SRC_FILES}
        ${IM_SOURCE_DIR}/randomscalarfield.cc
        ${IM_SOURCE_DIR}/randomvectorfield.cc
        ${IM_SOURCE_DIR}/gaussianscalar.cc
        ${IM_SOURCE_DIR}/lognormal.cc
        ${IM_SOURCE_DIR}/randomjf12.cc
        ${IM_SOURCE_DIR}/ensslinsteininger.cc
    )
endif()

add_library(tImagineModels SHARED ${IM_SRC_FILES})
target_include_directories(tImagineModels PUBLIC ${IM_INCLUDE_DIR}/ImagineModels)

if (FFTW_FOUND)
    target_include_directories(tImagineModels PUBLIC ${IM_INCLUDE_DIR}/ImagineModelsRandom)
endif()
target_link_libraries(tImagineModels PUBLIC ${LIBRARIES})
if (autodiff_FOUND) 
    target_include_directories(tImagineModels PUBLIC /usr/local/include/autodiff)
    target_include_directories(tImagineModels PUBLIC /usr/include/eigen3)
    target_link_libraries(tImagineModels PUBLIC autodiff::autodiff)
endif()

set_target_properties(tImagineModels PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${lib_path}" is_system)
if ("${is_system}" STREQUAL "-1")
    set_target_properties(tImagineModels PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
        INSTALL_RPATH "${lib_path}")
#    # The following is necessary for installation in a virtual
    # environment `python -m pip venv env`
    #s#et_target_properties(ImagineModels PROPERTIES
     #   INSTALL_RPATH_USE_LINK_PATH TRUE
     #   INSTALL_RPATH "${lib_path}")
        message("Hi internal: ")
endif()

enable_testing()
set(TESTSOURCES grid parameter_update positions)
foreach(test ${TESTSOURCES})
    set(test_name ${test}_ctest)
    add_executable(${test_name} "${PROJECT_SOURCE_DIR}/test/regular/test_${test}.cc")
    target_link_libraries(${test_name} tImagineModels ${LIBRARIES})
    add_test(${test_name} ${test_name})
endforeach()

configure_file(${PROJECT_SOURCE_DIR}/ImagineModels.pc.in ${CMAKE_BINARY_DIR}/ImagineModels.pc @ONLY)


message("CMAKE_INSTALL_INCLUDEDIR: "  ${CMAKE_INSTALL_INCLUDEDIR})
message("CMAKE_INSTALL_LIBDIR: "  ${CMAKE_INSTALL_LIBDIR})
message("CMAKE_INSTALL_RPATH: "  ${CMAKE_INSTALL_RPATH})
message("CMAKE_INSTALL_DATAROOTDIR: "  ${CMAKE_INSTALL_DATAROOTDIR})
install(TARGETS tImagineModels LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_BINARY_DIR}/ImagineModels.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)