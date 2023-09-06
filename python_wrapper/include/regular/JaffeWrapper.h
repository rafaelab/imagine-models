#ifndef JAFFEWRAPPER_H
#define JAFFEWRAPPER_H

#include <pybind11/pybind11.h>

#include "Jaffe.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Jaffe(py::module_ &m) {
    py::class_<JaffeParams>(m, "JaffeParams")
        .def(py::init<>())

        .def_readwrite("quadruple", &JaffeParams::quadruple)
        .def_readwrite("bss", &JaffeParams::bss)

        .def_readwrite("disk_amp", &JaffeParams::disk_amp)
        .def_readwrite("disk_z0", &JaffeParams::disk_z0)
        .def_readwrite("halo_amp",  &JaffeParams::halo_amp)
        .def_readwrite("halo_z0",  &JaffeParams::halo_z0)
        .def_readwrite("r_inner",  &JaffeParams::r_inner)
        .def_readwrite("r_scale",  &JaffeParams::r_scale)
        .def_readwrite("r_peak",  &JaffeParams::r_peak)

        .def_readwrite("ring", &JaffeParams::ring)
        .def_readwrite("bar", &JaffeParams::bar)
         // either ring or bar!
        .def_readwrite("ring_amp", &JaffeParams::ring_amp)
        .def_readwrite("ring_r", &JaffeParams::ring_r)
        .def_readwrite("bar_amp", &JaffeParams::bar_amp)
        .def_readwrite("bar_a", &JaffeParams::bar_a)
        .def_readwrite("bar_b", &JaffeParams::bar_b)
        .def_readwrite("bar_phi0", &JaffeParams::bar_phi0)

        .def_readwrite("arm_num",  &JaffeParams::arm_num)
        .def_readwrite("arm_r0",  &JaffeParams::arm_r0)
        .def_readwrite("arm_z0",  &JaffeParams::arm_z0)
        .def_readwrite("arm_phi1",  &JaffeParams::arm_phi1)
        .def_readwrite("arm_phi2",  &JaffeParams::arm_phi2)
        .def_readwrite("arm_phi3",  &JaffeParams::arm_phi3)
        .def_readwrite("arm_phi4",  &JaffeParams::arm_phi4)
        .def_readwrite("arm_amp1",  &JaffeParams::arm_amp1)
        .def_readwrite("arm_amp2",  &JaffeParams::arm_amp2)
        .def_readwrite("arm_amp3",  &JaffeParams::arm_amp3)
        .def_readwrite("arm_amp4",  &JaffeParams::arm_amp4)
        .def_readwrite("arm_pitch",  &JaffeParams::arm_pitch)

        .def_readwrite("comp_c",  &JaffeParams::comp_c)
        .def_readwrite("comp_d",  &JaffeParams::comp_d)
        .def_readwrite("comp_r",  &JaffeParams::comp_r)
        .def_readwrite("comp_p",  &JaffeParams::comp_p);

    py::class_<JaffeMagneticField, RegularVectorField>(m, "JaffeMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &JaffeMagneticField::at_position, "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
        .def_readwrite("param", &JaffeMagneticField::param);

            
}

#endif
