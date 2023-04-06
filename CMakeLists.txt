set(CMAKE_C_COMPILER "gcc-11")
set(CMAKE_CXX_COMPILER "/usr/bin/g++-11")

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")  # or 11, 14, 17, 20
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # optional, ensure standard is supported
set(CMAKE_CXX_EXTENSIONS ON)

add_compile_options(-std=c++17)

cmake_minimum_required(VERSION 3.4...3.18)
project(_ImagineModels)


include(GNUInstallDirs)

add_subdirectory(extern/pybind11)


pybind11_add_module(_ImagineModels  python_wrapper/ImagineModels.cc)

FIND_PACKAGE(autodiff)

message(autodiff_FOUND="${autodiff_FOUND}")

CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/c_library/headers/autodiff.hh.in ${CMAKE_SOURCE_DIR}/c_library/headers/autodiff.hh )


if (autodiff_FOUND)     
    set(autodiff_INCLUDE_DIRS /usr/include/eigen3)
    target_include_directories(_ImagineModels PUBLIC /usr/include/eigen3)
    find_package(Eigen3 REQUIRED NO_MODULE)
    message(Eigen3_FOUND="${Eigen3_FOUND}")

    get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
    foreach(dir ${dirs})
        message(STATUS "dir='${dir}'")
    endforeach()
    target_link_libraries(_ImagineModels PUBLIC autodiff::autodiff Eigen3::Eigen)


                
endif()

target_link_libraries(_ImagineModels PUBLIC "${CMAKE_SOURCE_DIR}/c_library/source/regularjf12.cc" 
                                             "${CMAKE_SOURCE_DIR}/c_library/source/jaffe.cc"
                                             "${CMAKE_SOURCE_DIR}/c_library/source/helix.cc"
                                             "${CMAKE_SOURCE_DIR}/c_library/source/ymw16.cc"
                                             "${CMAKE_SOURCE_DIR}/c_library/source/sun.cc"
                                             "${CMAKE_SOURCE_DIR}/c_library/source/tinyakovtkachev.cc"
                                             "${CMAKE_SOURCE_DIR}/c_library/source/hararimollerachroulet.cc"
                                             
                    )
                
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/c_library/findFFTW")

FIND_PACKAGE(FFTW)
message(FFTW_FOUND="${FFTW_FOUND}")
CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/c_library/headers/fftw.hh.in ${CMAKE_SOURCE_DIR}/c_library/headers/fftw.hh @ONLY)

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
target_compile_features(_ImagineModels PRIVATE cxx_std_17)
target_compile_definitions(_ImagineModels
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
