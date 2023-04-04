#include <pybind11/pybind11.h>

#include "../../../c_library/headers/Fauvet.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Fauvet(py::module_ &m) {
    py::class_<FauvetMagneticField, RegularVectorField>(m, "FauvetMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &FauvetMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("b_b0", &FauvetMagneticField::b_b0)
        .def_readwrite("b_r0", &FauvetMagneticField::b_r0)
        .def_readwrite("b_r_max", &FauvetMagneticField::b_r_max)
        .def_readwrite("b_r_min", &FauvetMagneticField::b_r_min)
        .def_readwrite("b_Rsun", &FauvetMagneticField::b_Rsun)
        .def_readwrite("b_z0", &FauvetMagneticField::b_z0)
        .def_readwrite("b_p", &FauvetMagneticField::b_p)
        .def_readwrite("b_chi0", &FauvetMagneticField::b_chi0)
        
        .def_readwrite("h_b0", &FauvetMagneticField::h_b0)
        .def_readwrite("h_r0", &FauvetMagneticField::h_r0)
        .def_readwrite("h_z0", &FauvetMagneticField::h_z0)
        .def_readwrite("h_z1a", &FauvetMagneticField::h_z1a)
        .def_readwrite("h_z1b", &FauvetMagneticField::h_z1b);
    
}
