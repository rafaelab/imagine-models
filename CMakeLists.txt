cmake_minimum_required(VERSION 3.4...3.18)
project(_ImagineModels)

add_subdirectory(extern/pybind11)

pybind11_add_module(_ImagineModels python_wrapper/imagine_models.cc)

target_link_libraries(_ImagineModels PUBLIC fftw3)

target_link_libraries(_ImagineModels PUBLIC "${CMAKE_SOURCE_DIR}/c_library/source/regularjf12.cc" 
                                            "${CMAKE_SOURCE_DIR}/c_library/source/jaffe.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/helix.cc"

                                            "${CMAKE_SOURCE_DIR}/c_library/source/randomjf12.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/ensslinsteininger.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/gaussianscalar.cc"
                                            "${CMAKE_SOURCE_DIR}/c_library/source/lognormal.cc"
# "${CMAKE_SOURCE_DIR}/c_library/source/ymw16.cc"
# "${CMAKE_SOURCE_DIR}/c_library/source/disc.cc"
# "${CMAKE_SOURCE_DIR}/c_library/source/component.cc"

 )


# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(_ImagineModels
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
