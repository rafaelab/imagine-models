#ifndef HELIXWRAPPER_H
#define HELIXWRAPPER_H

#include <pybind11/pybind11.h>

#include "Helix.h"

namespace py = pybind11;
using namespace pybind11::literals;


void Helix(py::module_ &m) {
    //py::implicitly_convertible<double, ad::real>()
    //py::implicitly_convertible<double, ad::detail::Real<1, double>>();
    
    py::class_<HelixParams>(m, "HelixParams")
        .def(py::init<>())
        .def_readwrite("ampx", &HelixParams::ampx)
        .def_readwrite("ampy", &HelixParams::ampy)
        .def_readwrite("ampz", &HelixParams::ampz)
        .def_readwrite("rmin", &HelixParams::rmin)
        .def_readwrite("rmax", &HelixParams::rmax);

    py::class_<HelixMagneticField, RegularVectorField>(m, "HelixMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", [](HelixMagneticField &self, double x, double y, double z)  {
            vector f = self.at_position(x, y, z);
            auto tp = std::make_tuple(f[0], f[1], f[2]);
            return tp;},
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership)
        .def("derivative", &HelixMagneticField::derivative,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("param", &HelixMagneticField::param);
}

#endif