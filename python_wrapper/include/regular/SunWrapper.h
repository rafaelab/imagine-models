#include <pybind11/pybind11.h>

#include "../../../c_library/headers/Sun.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Sun2008(py::module_ &m) {
    py::class_<Sun2008MagneticField, RegularVectorField>(m, "Sun2008MagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &Sun2008MagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("b_B0", &Sun2008MagneticField::b_B0)
        .def_readwrite("b_Rsun", &Sun2008MagneticField::b_Rsun)
        .def_readwrite("b_R0", &Sun2008MagneticField::b_R0)
        .def_readwrite("b_z0", &Sun2008MagneticField::b_z0)
        .def_readwrite("b_Rc", &Sun2008MagneticField::b_Rc)
        .def_readwrite("b_Bc", &Sun2008MagneticField::b_Bc)
        .def_readwrite("b_pitch_deg", &Sun2008MagneticField::b_pitch_deg)

        .def_readwrite("bH_B0", &Sun2008MagneticField::bH_B0)
        .def_readwrite("bH_z0", &Sun2008MagneticField::bH_z0)
        .def_readwrite("bH_z1a", &Sun2008MagneticField::bH_z1a)
        .def_readwrite("bH_z1b", &Sun2008MagneticField::bH_z1b)
        .def_readwrite("bH_R0", &Sun2008MagneticField::bH_R0);

}