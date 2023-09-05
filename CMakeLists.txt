cmake_minimum_required(VERSION 3.16)

project(_ImagineModels DESCRIPTION "A Python module for Galactic models.")

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


                
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/c_library/findFFTW")

                        )
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

include(CTest)
include(FetchContent)
include(GNUInstallDirs)
# find_package(PkgConfig REQUIRED)

# Activating the findFFTW cmake package
set(INTERNAL_findFFTW_DIR ${PROJECT_SOURCE_DIR}/extern/findFFTW)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${INTERNAL_findFFTW_DIR}")
if(NOT EXISTS ${INTERNAL_findFFTW_DIR})
    message(FATAL_ERROR "The findFFTW submodule was not initialized! Please run 'git submodule update --init'." )
endif()
find_package(FFTW QUIET COMPONENTS DOUBLE_LIB)

if(FFTW_FOUND)
    message("-- FFTW libraries found: ${FFTW_DOUBLE_LIB}; will include random field generation routines.")
        include_directories(${FFTW_INCLUDE_DIRS})
        set(LIBRARIES ${LIBRARIES} ${FFTW_DOUBLE_LIB})
else()
    message("-- FFTW libraries not found: random field generation disabled.")
endif()
add_compile_definitions(FFTW_FOUND=${FFTW_FOUND})
# pkg_search_module(autodiff QUIET)

add_subdirectory(${PROJECT_SOURCE_DIR}/c_library)

set(IM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/c_library/headers)

# Check for Python 3
set(INTERNAL_PYBIND11_DIR "${PROJECT_SOURCE_DIR}/extern/pybind11")
set(Python_FIND_VIRTUALENV FIRST)
find_package(Python3 COMPONENTS Interpreter Development)
if(Python3_FOUND)
    message("-- Python ${Python3_VERSION} found. Enabling support for pybind.")
    if(NOT EXISTS ${INTERNAL_PYBIND11_DIR})
        message(FATAL_ERROR "The pybind11 submodule was not initialized! Please run 'git submodule update --init'." )
    endif()
    add_subdirectory(${INTERNAL_PYBIND11_DIR})

    pybind11_add_module(_ImagineModels ${PROJECT_SOURCE_DIR}/python_wrapper/ImagineModels.cc)

    target_include_directories(_ImagineModels PUBLIC ${IM_INCLUDE_DIR}/ImagineModels)

    if(FFTW_FOUND)
        target_include_directories(_ImagineModels PUBLIC ${IM_INCLUDE_DIR}/ImagineModelsRandom)
        target_include_directories(_ImagineModels PUBLIC ${FFTW_INCLUDE_DIRS})
    endif()
    target_link_libraries(_ImagineModels PUBLIC ImagineModels ${LIBRARIES})
    # EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
    # define (VERSION_INFO) here.
    target_compile_definitions(_ImagineModels PRIVATE VERSION_INFO=${PROJECT_VERSION})
elseif(NOT Python3_FOUND)
    message("-- Could not find the Python Interpreter! Please provide provide a hint to CMake, using '-D CMAKE_PYTHON_EXECUTABLE'.")
else()
    message("-- Could not find the FFTW library. Disabling support for Python because of this.")
endif()