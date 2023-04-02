#ifndef FIELDBASES_H
#define FIELDBASES_H

#include <pybind11/pybind11.h>

#include "random_trampoline.h"
#include "regular_trampoline.h"

namespace py = pybind11;
using namespace pybind11::literals;

void FieldBases(py::module_ &m) {
    
    py::class_<Field<std::array<double, 3>, std::array<double*, 3>>,  PyVectorFieldBase>(m, "VectorFieldBase");

    py::class_<Field<double, double*>,  PyScalarFieldBase>(m, "ScalarFieldBase");

    py::class_<RandomField<std::array<double, 3>, std::array<double*, 3>>,  PyVectorRandomFieldBase>(m, "VectorRandomFieldBase");

    py::class_<RandomField<double, double*>,  PyScalarRandomFieldBase>(m, "ScalarRandomFieldBase");

}

#endif