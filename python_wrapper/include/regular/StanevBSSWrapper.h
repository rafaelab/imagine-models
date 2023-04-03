#include <pybind11/pybind11.h>

#include "../../../c_library/headers/StanevBSS.h"

namespace py = pybind11;
using namespace pybind11::literals;

void StanevBSS(py::module_ &m) {
    py::class_<StanevBSSMagneticField, RegularVectorField>(m, "StanevBSSMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &StanevBSSMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("b_b0", &StanevBSSMagneticField::b_b0)
        .def_readwrite("b_r0", &StanevBSSMagneticField::b_r0)
        .def_readwrite("b_phi0", &StanevBSSMagneticField::b_phi0)
        .def_readwrite("b_r_max", &StanevBSSMagneticField::b_r_max)
        .def_readwrite("b_r_min", &StanevBSSMagneticField::b_r_min)
        .def_readwrite("b_Rsun", &StanevBSSMagneticField::b_Rsun)
        .def_readwrite("b_z0", &StanevBSSMagneticField::b_z0)
        .def_readwrite("b_p", &StanevBSSMagneticField::b_p);
    
}