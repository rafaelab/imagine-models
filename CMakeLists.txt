cmake_minimum_required(VERSION 3.4...3.18)
project(ImagineModels)

add_subdirectory(extern/pybind11)

pybind11_add_module(ImagineModels ImagineModels/python/imagine_models.cc)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(ImagineModels
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
