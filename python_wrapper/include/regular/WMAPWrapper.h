#include <pybind11/pybind11.h>

#include "../../../c_library/headers/WMAP.h"

namespace py = pybind11;
using namespace pybind11::literals;

void WMAP(py::module_ &m) {
    py::class_<WMAPMagneticField, RegularVectorField>(m, "WMAPMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &WMAPMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("b_b0", &WMAPMagneticField::b_b0)
        .def_readwrite("b_r0", &WMAPMagneticField::b_r0)
        .def_readwrite("b_r_max", &WMAPMagneticField::b_r_max)
        .def_readwrite("b_r_min", &WMAPMagneticField::b_r_min)
        .def_readwrite("b_Rsun", &WMAPMagneticField::b_Rsun)
        .def_readwrite("b_z0", &WMAPMagneticField::b_z0)
        
        .def_readwrite("b_psi0", &WMAPMagneticField::b_psi0)
        .def_readwrite("b_psi1", &WMAPMagneticField::b_psi1)
        .def_readwrite("b_xsi0", &WMAPMagneticField::b_xsi0)
        .def_readwrite("b_anti", &WMAPMagneticField::anti);
    
}