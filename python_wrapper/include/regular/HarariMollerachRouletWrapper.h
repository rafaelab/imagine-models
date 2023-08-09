#include <pybind11/pybind11.h>

#include "HarariMollerachRoulet.h"

namespace py = pybind11;
using namespace pybind11::literals;

void HarariMollerachRoulet(py::module_ &m) {
    py::class_<HMRMagneticField, RegularVectorField>(m, "HMRMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &HMRMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("b_Rsun", &HMRMagneticField::b_Rsun)
        .def_readwrite("b_r_max", &HMRMagneticField::b_r_max)
        .def_readwrite("b_z1", &HMRMagneticField::b_z1)
        .def_readwrite("b_z2", &HMRMagneticField::b_z2)
        .def_readwrite("b_r1", &HMRMagneticField::b_r1)
        .def_readwrite("b_p", &HMRMagneticField::b_p)
        .def_readwrite("b_epsilon0", &HMRMagneticField::b_epsilon0);
}