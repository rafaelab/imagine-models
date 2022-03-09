#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "trampoline.h"

namespace py = pybind11;
using namespace pybind11::literals;

// Trampoline classes via templates



PYBIND11_MODULE(ImagineModels, m) {
    m.doc() = "IMAGINE Magnetic Field Model Library";

    py::class_<MagneticField, PyMagneticField>(m, "Magnetic_Field")
        .def(py::init<>());

    py::class_<RegularMagneticField, MagneticField, PyRegularMagneticField>(m, "RegularMagneticField")
        .def(py::init<>())
        .def("evaluate_at_pos", &RegularMagneticField::evaluate_at_pos, "pos"_a)
        .def("evaluate_grid", &RegularMagneticField::evaluate_at_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a);

   py::class_<JF12MagneticField, RegularMagneticField, PyJF12MagneticField>(m, "JF12MagneticField")
        .def(py::init<>())
        .def("evaluate_at_pos", &JF12MagneticField::evaluate_at_pos, "pos"_a);
    }
