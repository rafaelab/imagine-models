#include <pybind11/pybind11.h>
#include <autodiff/forward/real.hpp>

#include "../../../c_library/headers/Helix.h"
#include "../autodiff_wrapper.h"

namespace py = pybind11;
using namespace pybind11::literals;
namespace ad = autodiff;




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

        .def("at_position", &HelixMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
        .def("derivative", &HelixMagneticField::derivative,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("param", &HelixMagneticField::param);
}
