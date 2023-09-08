#ifndef UNIFORMWRAPPER_H
#define UNIFORMWRAPPER_H

#include <pybind11/pybind11.h>

#include "Uniform.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Uniform(py::module_ &m) {

py::class_<UniformParams>(m, "UniformParams")
        .def(py::init<>())
        .def_readwrite("bx", &UniformParams::bx)
        .def_readwrite("by", &UniformParams::by)
        .def_readwrite("bz", &UniformParams::bz);


py::class_<UniformMagneticField, RegularVectorField>(m, "UniformMagneticField")
    .def(py::init<>())
    .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
    .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

    .def("at_position", [](UniformMagneticField &self, double x, double y, double z)  {
        vector f = self.at_position(x, y, z);
        auto tp = std::make_tuple(f[0], f[1], f[2]);
        return tp;},
        "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership)

    .def_readwrite("param", &UniformMagneticField::param);
}

#endif
