project(ImagineModels VERSION 0.1 DESCRIPTION "A library for Galactic models.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

set(IM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/headers/)
include_directories(${IM_INCLUDE_DIR})
set(IM_SRC_DIR ${PROJECT_SOURCE_DIR}/source/)
set(IM_SRC_FILES
    ${IM_SRC_DIR}/helix.cc
    ${IM_SRC_DIR}/jaffe.cc
    ${IM_SRC_DIR}/regularjf12.cc
    ${IM_SRC_DIR}/ymw16.cc
    ${IM_SRC_DIR}/sun.cc
    ${IM_SRC_DIR}/wmap.cc
    ${IM_SRC_DIR}/stanevbss.cc
    ${IM_SRC_DIR}/fauvet.cc
    ${IM_SRC_DIR}/tinyakovtkachev.cc
    ${IM_SRC_DIR}/hararimollerachroulet.cc   
)

if(FFTW_FOUND)
    set(IM_SRC_FILES
        ${IM_SRC_FILES}
        ${IM_SRC_DIR}/gaussianscalar.cc
        ${IM_SRC_DIR}/lognormal.cc
        ${IM_SRC_DIR}/randomjf12.cc
        ${IM_SRC_DIR}/ensslinsteininger.cc
    )
endif()

add_library(ImagineModels SHARED ${IM_SRC_FILES})
target_link_libraries(ImagineModels PUBLIC ${LIBRARIES})
set_target_properties(ImagineModels PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

enable_testing()
set(TESTSOURCES grid parameter_update positions)
foreach(test ${TESTSOURCES})
    set(test_name ${test}_ctest)
    add_executable(${test_name} "${PROJECT_SOURCE_DIR}/test/regular/test_${test}.cc")
    target_link_libraries(${test_name} ImagineModels ${LIBRARIES})
    add_test(${test_name} ${test_name})
endforeach()

configure_file(ImagineModels.pc.in ImagineModels.pc @ONLY)

install(TARGETS ImagineModels LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_BINARY_DIR}/ImagineModels.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)