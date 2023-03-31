#include <pybind11/pybind11.h>

#include "../../../c_library/headers/Helix.h"

namespace py = pybind11;
using namespace pybind11::literals;


void Helix(py::module_ &m) {
    py::class_<HelixMagneticField, RegularVectorField>(m, "HelixMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &HelixMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("ampx", &HelixMagneticField::ampx)
        .def_readwrite("ampy", &HelixMagneticField::ampy)
        .def_readwrite("ampz", &HelixMagneticField::ampz)
        .def_readwrite("rmin", &HelixMagneticField::rmin)
        .def_readwrite("rmax", &HelixMagneticField::rmax);
}
