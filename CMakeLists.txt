cmake_minimum_required(VERSION 3.4...3.18)
project(ImagineModels)

add_subdirectory(extern/pybind11)

pybind11_add_module(ImagineModels ImagineModels/python_wrapper/imagine_models.cc)
target_link_libraries(ImagineModels PUBLIC "${CMAKE_SOURCE_DIR}/ImagineModels/c/source/jf12.cc"
                                           "${CMAKE_SOURCE_DIR}/ImagineModels/c/source/regular.cc")

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(ImagineModels
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
