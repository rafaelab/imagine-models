#include "../../../c_library/headers/RegularJF12.h"


void RegularJF12(py::module_ &m) {
    py::class_<JF12RegularParams>(m, "JF12RegularParams")
        .def(py::init<>())
        .def_readwrite("b_arm_1", &JF12RegularParams::b_arm_1)
        .def_readwrite("b_arm_2", &JF12RegularParams::b_arm_2)
        .def_readwrite("b_arm_3", &JF12RegularParams::b_arm_3)
        .def_readwrite("b_arm_4", &JF12RegularParams::b_arm_4)
        .def_readwrite("b_arm_5", &JF12RegularParams::b_arm_5)
        .def_readwrite("b_arm_6", &JF12RegularParams::b_arm_6)
        .def_readwrite("b_arm_7", &JF12RegularParams::b_arm_7)
        .def_readwrite("b_ring", &JF12RegularParams::b_ring)
        .def_readwrite("h_disk", &JF12RegularParams::h_disk)
        .def_readwrite("w_disk", &JF12RegularParams::w_disk)
        .def_readwrite("Bn", &JF12RegularParams::Bn)
        .def_readwrite("Bs", &JF12RegularParams::Bs)
        .def_readwrite("rn", &JF12RegularParams::rn)
        .def_readwrite("rs", &JF12RegularParams::rs)
        .def_readwrite("wh", &JF12RegularParams::wh)
        .def_readwrite("z0", &JF12RegularParams::z0)
        .def_readwrite("B0_X", &JF12RegularParams::B0_X)
        .def_readwrite("Xtheta_const", &JF12RegularParams::Xtheta_const)
        .def_readwrite("rpc_X", &JF12RegularParams::rpc_X)
        .def_readwrite("r0_X", &JF12RegularParams::r0_X);


    py::class_<JF12MagneticField, RegularVectorField>(m, "JF12RegularField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &JF12MagneticField::at_position, "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
        .def("derivative", &JF12MagneticField::derivative,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("param", &JF12MagneticField::param);

}