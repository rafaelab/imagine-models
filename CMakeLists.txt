cmake_minimum_required(VERSION 3.4...3.18)
project(_ImagineModels)

CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/c_library/include/config/autodiff.hh.in ${CMAKE_SOURCE_DIR}/c_library/include/config/autodiff.hh )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/c_library/findFFTW")


FIND_PACKAGE(FFTW)

message(CMAKE_MODULE_PATH="${CMAKE_MODULE_PATH}")
message(PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
message(FFTW_FOUND="${FFTW_FOUND}")


CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/c_library/include/config/fftw.hh.in ${CMAKE_SOURCE_DIR}/c_library/include/config/fftw.hh @ONLY)


add_subdirectory(extern/pybind11)

pybind11_add_module(_ImagineModels python_wrapper/ImagineModels.cc)


target_link_libraries(_ImagineModels PUBLIC "${CMAKE_SOURCE_DIR}/c_library/source/regularjf12.cc" 
                                            "${CMAKE_SOURCE_DIR}/c_library/source/jaffe.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/helix.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/ymw16.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/sun.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/tinyakovtkachev.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/hararimollerachroulet.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/stanevbss.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/wmap.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/fauvet.cc"
                                             
                    )

if (FFTW_FOUND) 
    target_link_libraries(_ImagineModels PUBLIC fftw3)

    target_link_libraries(_ImagineModels PUBLIC "${CMAKE_SOURCE_DIR}/c_library/source/randomjf12.cc"
                                                "${CMAKE_SOURCE_DIR}/c_library/source/ensslinsteininger.cc"
                                                "${CMAKE_SOURCE_DIR}/c_library/source/gaussianscalar.cc"
                                                "${CMAKE_SOURCE_DIR}/c_library/source/lognormal.cc"
                        )
endif()



# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(_ImagineModels
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
