#include <pybind11/pybind11.h>

#include "GaussianScalar.h"

namespace py = pybind11;
using namespace pybind11::literals;

void GaussianScalar(py::module_ &m) {
    py::class_<GaussianScalarField, RandomScalarField>(m, "GaussianScalarField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("mean", &GaussianScalarField::mean)
        .def_readwrite("spectral_amplitude", &GaussianScalarField::spectral_amplitude)
        .def_readwrite("spectral_offset", &GaussianScalarField::spectral_offset)
        .def_readwrite("spectral_slope", &GaussianScalarField::spectral_slope);
}