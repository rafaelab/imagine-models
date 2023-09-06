#ifndef HMRWRAPPER_H
#define HMRWRAPPER_H

#include <pybind11/pybind11.h>

#include "HarariMollerachRoulet.h"

namespace py = pybind11;
using namespace pybind11::literals;

void HarariMollerachRoulet(py::module_ &m) {
    py::class_<HMRParams>(m, "HMRParams")
        .def(py::init<>())
        .def_readwrite("b_Rsun", &HMRParams::b_Rsun)
        .def_readwrite("b_r_max", &HMRParams::b_r_max)
        .def_readwrite("b_z1", &HMRParams::b_z1)
        .def_readwrite("b_z2", &HMRParams::b_z2)
        .def_readwrite("b_r1", &HMRParams::b_r1)
        .def_readwrite("b_p", &HMRParams::b_p)
        .def_readwrite("b_epsilon0", &HMRParams::b_epsilon0);

    py::class_<HMRMagneticField, RegularVectorField>(m, "HMRMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &HMRMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("param", &HMRMagneticField::param);


}

#endif