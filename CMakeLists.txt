cmake_minimum_required(VERSION 3.16)

project(_ImagineModels DESCRIPTION "A Python module for Galactic models.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

include(CTest)
include(FetchContent)
include(GNUInstallDirs)


# The follwing line of code are form 
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
#set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#message("CMAKE_INSTALL_RPATH: "  ${CMAKE_INSTALL_RPATH})
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

## end RPATH code

set(USE_AUTODIFF ON CACHE BOOL OFF)
set(USE_FFTW ON CACHE BOOL OFF)

if (USE_AUTODIFF)
    message("-- Trying to find autodiff. Turn off via 'cmake ..DUSE_AUTODIFF=OFF'")
    find_package(autodiff)
    if (autodiff_FOUND)         
        message("-- autodiff library found: attempting to find Eigen3.")
        set(autodiff_INCLUDE_DIRS /usr/include/eigen3)
    
        find_package(Eigen3 REQUIRED NO_MODULE)

        if(NOT Eigen3_FOUND)
            message("-- Could not find the Eigen3 library! Autodifferention disabled.")
        else()
            message("-- Eigen3 library found: Autodifferention enabled.")
    endif()
    else()
        message("-- Could not find the autodiff library! Autodifferention disabled.")
    endif()
else()
    set(autodiff_FOUND OFF) 
    message("-- autodiff library manually disabled")
endif()

add_compile_definitions(autodiff_FOUND=${autodiff_FOUND})
                
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/c_library/findFFTW")

# find_package(PkgConfig REQUIRED)

# Activating the findFFTW cmake package
if (USE_FFTW)
    message("-- Trying to find FFTW. Turn off via 'cmake ..DUSE_FFTW=OFF'")
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
else()
    set(FFTW_FOUND OFF) 
    message("-- fftw library manually disabled")
endif()

if (FFTW_FOUND)
    add_compile_definitions(FFTW_FOUND=1)
else()
    add_compile_definitions(FFTW_FOUND=0)
endif()
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
    if (autodiff_FOUND) 
        target_include_directories(_ImagineModels PUBLIC /usr/local/include/autodiff)
        target_include_directories(_ImagineModels PUBLIC /usr/include/eigen3)
    endif()

    set_target_properties(_ImagineModels PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

    # target_link_libraries(_ImagineModels PRIVATE supporting_library)
    
    set_target_properties(_ImagineModels PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

    message("LIBRARIES: "  ${LIBRARIES})
    target_link_libraries(_ImagineModels PUBLIC tImagineModels ${LIBRARIES})
    if (autodiff_FOUND)
        target_link_libraries(_ImagineModels PUBLIC autodiff::autodiff)
    endif()

    
   
    if(UNIX AND NOT APPLE)
        message("unix: ${Python3_SITEARCH}") 
        set_target_properties(_ImagineModels PROPERTIES INSTALL_RPATH "${Python3_SITEARCH}")
    elseif(APPLE)
        message("apple: @loader_path/site_packages")
        set_target_properties(_ImagineModels PROPERTIES INSTALL_RPATH "@loader_path/site_packages")
    else() 
        message("No unix and no apple")
    endif()
    
    # EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
    # define (VERSION_INFO) here.
    target_compile_definitions(_ImagineModels PRIVATE VERSION_INFO=${PROJECT_VERSION})
    set(install_path "${Python3_SITEARCH}")
    # string(REPLACE "python3" "" pypath ppath)
    message("Python exec ${Python3_EXECUTABLE}")
    # set(lib_path "${Python3_EXECUTABLE}/${CMAKE_INSTALL_LIBDIR}")
    # message( "Python3_pypath: ${ppath}")
    message( "Python3_STDARCH: ${Python3_STDARCH}")
    message( "Python3_SITEARCH: ${Python3_SITEARCH}")
    message( "Python3_BINARCH: ${Python3_BINARCH}")
    message( "Python3_LIBRARIES: ${Python3_LIBRARIES}")
    message( "Python3_LIBRARY_DIRS: ${Python3_LIBRARY_DIRS}")
    message( "Python3_RUNTIME_LIBRARY_DIRS: ${Python3_RUNTIME_LIBRARY_DIRS}")
    
else()
    message("-- Could not find the Python Interpreter! If you want to build the python package, please provide provide a hint to CMake, using '-D CMAKE_PYTHON_EXECUTABLE'.")
endif()

#if (SKBUILD)
#    find_package(PythonExtensions REQUIRED)
#    set(lib_path "${PYTHON_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
#else()
# set(lib_path "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
#endif()
#list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${install_path}" is_system)
#if ("${is_system}" STREQUAL "-1")
#    set_target_properties(_ImagineModels PROPERTIES
#        INSTALL_RPATH_USE_LINK_PATH TRUE
#        INSTALL_RPATH "${install_path}")
#    # The following is necessary for installation in a virtual
#    # environment `python -m pip venv env`
#    set_target_properties(_ImagineModels PROPERTIES
#        INSTALL_RPATH_USE_LINK_PATH TRUE
#        INSTALL_RPATH "${install_path}")
#        message("Hi: ")
#endif()

#message("install_path: "  ${install_path})
#message("is_system: "  ${is_system})