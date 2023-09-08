#ifndef TTWRAPPER_H
#define TTWRAPPER_H

#include <pybind11/pybind11.h>

#include "TinyakovTkachev.h"

namespace py = pybind11;
using namespace pybind11::literals;

void TinyakovTkachev(py::module_ &m) {
    py::class_<TTParams>(m, "TTParams")
        .def(py::init<>())
        .def_readwrite("b_b0", &TTParams::b_b0)
        .def_readwrite("b_Rsun", &TTParams::b_Rsun)
        .def_readwrite("b_r_max", &TTParams::b_r_max)
        .def_readwrite("b_r_min", &TTParams::b_r_min)
        .def_readwrite("b_d", &TTParams::b_d)
        .def_readwrite("b_z0", &TTParams::b_z0)
        .def_readwrite("b_p", &TTParams::b_p);

    py::class_<TTMagneticField, RegularVectorField>(m, "TTMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", [](TTMagneticField &self, double x, double y, double z)  {
            vector f = self.at_position(x, y, z);
            auto tp = std::make_tuple(f[0], f[1], f[2]);
            return tp;},
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership)

        .def_readwrite("param", &TTMagneticField::param);

}

#endif