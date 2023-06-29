cmake_minimum_required(VERSION 3.18)

project(_ImagineModels DESCRIPTION "A Python module for Galactic models.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

include(CTest)
include(FetchContent)
include(GNUInstallDirs)
# find_package(PkgConfig REQUIRED)

# Fetching the findFFTW cmake package
set(INTERNAL_findFFTW_DIR ${PROJECT_SOURCE_DIR}/extern/findFFTW)
FetchContent_Declare(findFFTW
    GIT_REPOSITORY https://github.com/egpbos/findFFTW
    GIT_TAG ac788529392825067c9118f9f099aaa1efb77589
    GIT_SHALLOW ON
    GIT_SUBMODULES ""
    GIT_CONFIG advice.detachedHead=False
    SOURCE_DIR ${INTERNAL_findFFTW_DIR}
)
FetchContent_Populate(findFFTW)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${INTERNAL_findFFTW_DIR}")

find_package(FFTW QUIET COMPONENTS DOUBLE_LIB)
if(FFTW_FOUND)
message("-- FFTW libraries found: ${FFTW_DOUBLE_LIB}; will include random field generation routines.")
    include_directories(${FFTW_INCLUDE_DIRS})
    set(LIBRARIES ${LIBRARIES} ${FFTW_DOUBLE_LIB})
endif()
add_compile_definitions(FFTW_FOUND="${FFTW_FOUND}")
# pkg_search_module(autodiff QUIET)

add_subdirectory(${PROJECT_SOURCE_DIR}/c_library)

set(IM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/c_library/headers/)
include_directories(${IM_INCLUDE_DIR})

# Check for Python 3
set(INTERNAL_PYBIND11_DIR "${PROJECT_SOURCE_DIR}/extern/pybind11")
set(Python_FIND_VIRTUALENV FIRST)
find_package(Python3 COMPONENTS Interpreter Development)
if(Python3_FOUND AND FFTW_FOUND)
    message("-- Python ${Python3_VERSION} and FFTW ${FFTW_VERSION} found. Enabling support for pybind.")
    if(NOT EXISTS ${INTERNAL_PYBIND11_DIR})
        message(FATAL_ERROR "Pybind11 submodule was not materialized!" )
    endif()
    # set(pybind11_INCLUDE_DIR ${INTERNAL_PYBIND11_DIR}/include)
    add_subdirectory(${INTERNAL_PYBIND11_DIR})
    pybind11_add_module(_ImagineModels ${PROJECT_SOURCE_DIR}/python_wrapper/ImagineModels.cc)
    target_include_directories(_ImagineModels PUBLIC ${FFTW_INCLUDE_DIRS})
    target_link_libraries(_ImagineModels PUBLIC ImagineModels ${LIBRARIES})
    # EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
    # define (VERSION_INFO) here.
    target_compile_definitions(_ImagineModels PRIVATE VERSION_INFO=${PROJECT_VERSION})
elseif(NOT Python3_FOUND)
    message("-- Could not find the Python Interpreter! Please provide provide a hint to CMake, using '-D CMAKE_PYTHON_EXECUTABLE'.")
else()
    message("-- Could not find the FFTW library. Disabling support for Python because of this.")
endif()