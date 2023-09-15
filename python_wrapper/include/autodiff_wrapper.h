#include <iostream>
#include <pybind11/pybind11.h>
#include <autodiff/forward/real.hpp>


namespace py = pybind11;
namespace ad = autodiff;


namespace pybind11 { namespace detail {

template <> struct type_caster<ad::real> : public type_caster_base<ad::real> {
    using base = type_caster_base<ad::real>;
public:
    bool load(handle src, bool convert) {
        if (py::isinstance<py::float_>(src)) {
            value = new ad::real(py::cast<double>(src));
            return true;
        }
        return false;
    }

    static handle cast(ad::real src, return_value_policy policy, handle parent) {
        return PyFloat_FromDouble(src.val());
    }
};

}}