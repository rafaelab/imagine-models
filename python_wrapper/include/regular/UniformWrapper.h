#include <pybind11/pybind11.h>

#include "Uniform.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Uniform(py::module_ &m) {
py::class_<UniformMagneticField, RegularVectorField>(m, "UniformMagneticField")
    .def(py::init<>())
    .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
    .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

    .def("at_position", &UniformMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

    .def_readwrite("bx", &UniformMagneticField::bx)
    .def_readwrite("by", &UniformMagneticField::by)
    .def_readwrite("bz", &UniformMagneticField::bz);
}
