#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

#include "trampoline.h"

namespace py = pybind11;
using namespace pybind11::literals;

// Trampoline classes via templates


PYBIND11_MODULE(_ImagineModels, m) {
    m.doc() = "IMAGINE Magnetic Field Model Library";

    py::class_<MagneticField, PyMagneticField>(m, "MagneticField")
        .def(py::init<>());

    py::class_<RegularMagneticField, MagneticField, PyRegularMagneticField>(m, "RegularMagneticField")
        .def(py::init<>())
        .def("_evaluate_grid", &RegularMagneticField::_evaluate_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a, "ev_at_pos"_a)
        .def("evaluate_grid", &RegularMagneticField::evaluate_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a);

   py::class_<JF12MagneticField, RegularMagneticField, PyJF12MagneticField>(m, "JF12MagneticField")
        .def(py::init<>())
        .def("evaluate_at_pos", &JF12MagneticField::evaluate_at_pos, "pos"_a)
        .def_readwrite("b_arm_1", &JF12MagneticField::b_arm_1)
        .def_readwrite("b_arm_2", &JF12MagneticField::b_arm_2)
        .def_readwrite("b_arm_3", &JF12MagneticField::b_arm_3)
        .def_readwrite("b_arm_4", &JF12MagneticField::b_arm_4)
        .def_readwrite("b_arm_5", &JF12MagneticField::b_arm_5)
        .def_readwrite("b_arm_6", &JF12MagneticField::b_arm_6)
        .def_readwrite("b_arm_7", &JF12MagneticField::b_arm_7)
        .def_readwrite("b_ring", &JF12MagneticField::b_ring)
        .def_readwrite("h_disk", &JF12MagneticField::h_disk)
        .def_readwrite("w_disk", &JF12MagneticField::w_disk)
        .def_readwrite("Bn", &JF12MagneticField::Bn)
        .def_readwrite("Bs", &JF12MagneticField::Bs)
        .def_readwrite("rn", &JF12MagneticField::rn)
        .def_readwrite("rs", &JF12MagneticField::rs)
        .def_readwrite("wh", &JF12MagneticField::wh)
        .def_readwrite("z0", &JF12MagneticField::z0)
        .def_readwrite("B0_X", &JF12MagneticField::B0_X)
        .def_readwrite("Xtheta_const", &JF12MagneticField::Xtheta_const)
        .def_readwrite("rpc_X", &JF12MagneticField::rpc_X)
        .def_readwrite("r0_X", &JF12MagneticField::r0_X);
    }
