#ifndef SUNWRAPPER_H
#define SUNWRAPPER_H

#include <pybind11/pybind11.h>

#include "Sun.h"

void Sun2008(py::module_ &m) {

    py::class_<SunParams>(m, "SunParams")
        .def(py::init<>())
        .def_readwrite("b_B0", &SunParams::b_B0)
        .def_readwrite("b_Rsun", &SunParams::b_Rsun)
        .def_readwrite("b_R0", &SunParams::b_R0)
        .def_readwrite("b_z0", &SunParams::b_z0)
        .def_readwrite("b_Rc", &SunParams::b_Rc)
        .def_readwrite("b_Bc", &SunParams::b_Bc)
        .def_readwrite("b_pitch_deg", &SunParams::b_pitch_deg)

        .def_readwrite("bH_B0", &SunParams::bH_B0)
        .def_readwrite("bH_z0", &SunParams::bH_z0)
        .def_readwrite("bH_z1a", &SunParams::bH_z1a)
        .def_readwrite("bH_z1b", &SunParams::bH_z1b)
        .def_readwrite("bH_R0", &SunParams::bH_R0);

    py::class_<SunMagneticField, RegularVectorField>(m, "SunMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &SunMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("param", &SunMagneticField::param);
}

#endif