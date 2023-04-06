#include "../../../c_library/headers/TinyakovTkachev.h"

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

        .def("at_position", &TTMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("param", &TTMagneticField::param);

}