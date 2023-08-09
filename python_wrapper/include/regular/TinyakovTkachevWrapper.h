#include <pybind11/pybind11.h>

#include "ImagineModels/TinyakovTkachev.h"

namespace py = pybind11;
using namespace pybind11::literals;

void TinyakovTkachev(py::module_ &m) {
    py::class_<TTMagneticField, RegularVectorField>(m, "TTMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &TTMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("b_b0", &TTMagneticField::b_b0)
        .def_readwrite("b_Rsun", &TTMagneticField::b_Rsun)
        .def_readwrite("b_r_max", &TTMagneticField::b_r_max)
        .def_readwrite("b_r_min", &TTMagneticField::b_r_min)
        .def_readwrite("b_d", &TTMagneticField::b_d)
        .def_readwrite("b_z0", &TTMagneticField::b_z0)
        .def_readwrite("b_p", &TTMagneticField::b_p);
    
}