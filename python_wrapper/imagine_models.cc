#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

#include "magnetic_trampoline.h"
#include "thermal_trampoline.h"

namespace py = pybind11;
using namespace pybind11::literals;


PYBIND11_MODULE(_ImagineModels, m) {
    m.doc() = "IMAGINE Model Library";

/////////////////////////////////Magnetic Field/////////////////////////////////

// Base Class
    py::class_<MagneticField, PyMagneticField>(m, "MagneticField")
        .def(py::init<>());

//Derived Base class for Regular field

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

    py::class_<HelixMagneticField, RegularMagneticField, PyHelixMagneticField>(m, "HelixMagneticField")
        .def(py::init<>())
        .def("evaluate_at_pos", &HelixMagneticField::evaluate_at_pos, "pos"_a)
        .def("dampx_grid", &HelixMagneticField::dampx_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a)
        .def("dampy_grid", &HelixMagneticField::dampy_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a)
        .def("dampz_grid", &HelixMagneticField::dampz_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a)
        .def_readwrite("ampx", &HelixMagneticField::ampx)
        .def_readwrite("ampy", &HelixMagneticField::ampy)
        .def_readwrite("ampz", &HelixMagneticField::ampz)
        .def_readwrite("rmin", &HelixMagneticField::rmin)
        .def_readwrite("rmax", &HelixMagneticField::rmax);

    py::class_<JaffeMagneticField, RegularMagneticField, PyJaffeMagneticField>(m, "JaffeMagneticField")
        .def(py::init<>())
        .def("evaluate_at_pos", &JaffeMagneticField::evaluate_at_pos, "pos"_a)
        .def_readwrite("quadruple", &JaffeMagneticField::quadruple)
        .def_readwrite("bss", &JaffeMagneticField::bss)

        .def_readwrite("disk_amp", &JaffeMagneticField::disk_amp)
        .def_readwrite("disk_z0", &JaffeMagneticField::disk_z0)
        .def_readwrite("halo_amp",  &JaffeMagneticField::halo_amp)
        .def_readwrite("halo_z0",  &JaffeMagneticField::halo_z0)
        .def_readwrite("r_inner",  &JaffeMagneticField::r_inner)
        .def_readwrite("r_scale",  &JaffeMagneticField::r_scale)
        .def_readwrite("r_peak",  &JaffeMagneticField::r_peak)

        .def_readwrite("ring", &JaffeMagneticField::ring)
        .def_readwrite("bar", &JaffeMagneticField::bar)
         // either ring or bar!
        .def_readwrite("ring_amp", &JaffeMagneticField::ring_amp)
        .def_readwrite("ring_r", &JaffeMagneticField::ring_r)
        .def_readwrite("bar_amp", &JaffeMagneticField::bar_amp)
        .def_readwrite("bar_a", &JaffeMagneticField::bar_a)
        .def_readwrite("bar_b", &JaffeMagneticField::bar_b)
        .def_readwrite("bar_phi0", &JaffeMagneticField::bar_phi0)

        .def_readwrite("arm_num",  &JaffeMagneticField::arm_num)
        .def_readwrite("arm_r0",  &JaffeMagneticField::arm_r0)
        .def_readwrite("arm_z0",  &JaffeMagneticField::arm_z0)
        .def_readwrite("arm_phi1",  &JaffeMagneticField::arm_phi1)
        .def_readwrite("arm_phi2",  &JaffeMagneticField::arm_phi2)
        .def_readwrite("arm_phi3",  &JaffeMagneticField::arm_phi3)
        .def_readwrite("arm_phi4",  &JaffeMagneticField::arm_phi4)
        .def_readwrite("arm_amp1",  &JaffeMagneticField::arm_amp1)
        .def_readwrite("arm_amp2",  &JaffeMagneticField::arm_amp2)
        .def_readwrite("arm_amp3",  &JaffeMagneticField::arm_amp3)
        .def_readwrite("arm_amp4",  &JaffeMagneticField::arm_amp4)
        .def_readwrite("arm_pitch",  &JaffeMagneticField::arm_pitch)

        .def_readwrite("comp_c",  &JaffeMagneticField::comp_c)
        .def_readwrite("comp_d",  &JaffeMagneticField::comp_d)
        .def_readwrite("comp_r",  &JaffeMagneticField::comp_r)
        .def_readwrite("comp_p",  &JaffeMagneticField::comp_p);

/////////////////////////////Thermal Electron Field/////////////////////////////

    py::class_<ThermalElectronField, PyThermalElectronField>(m, "ThermalElectronField")
        .def(py::init<>());
    }

    py::class_<RegularThermalElectronField, ThermalElectronField, PyRegularThermalElectronField>(m, "RegularThermalElectronField")
        .def(py::init<>())
        .def("_evaluate_grid", &RegularThermalElectronField::_evaluate_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a, "ev_at_pos"_a)
        .def("evaluate_grid", &RegularThermalElectronField::evaluate_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a);

    py::class_<YMW16ThickDisc, YMW16Component, PyYMW16ThickDisc>(m, "YMW16ThickDisc")
        .def(py::init<>())
        .def("evaluate_at_pos", &YMW16ThickDisc::evaluate_at_pos, "pos"_a)

        .def_readwrite("t1_ad",  &YMW16ThickDisc::t1_ad)
        .def_readwrite("t1_bd",  &YMW16ThickDisc::t1_bd)
        .def_readwrite("t1_n1",  &YMW16ThickDisc::t1_n1)
        .def_readwrite("t1_h1",  &YMW16ThickDisc::t1_h1);

   py::class_<YMW16ThinDisc, YMW16Component, PyYMW16ThinDisc>(m, "YMW16ThinDisc")
        .def(py::init<>())
        .def("evaluate_at_pos", &YMW16ThinDisc::evaluate_at_pos, "pos"_a)

        .def_readwrite("t1_ad",  &YMW16ThinDisc::t1_ad)
        .def_readwrite("t1_bd",  &YMW16ThinDisc::t1_bd)
        .def_readwrite("t2_n2",  &YMW16ThinDisc::t2_n2)
        .def_readwrite("t2_k2",  &YMW16ThinDisc::t2_h2);
        .def_readwrite("t2_a2",  &YMW16ThinDisc::t2_a2)
        .def_readwrite("t2_b2",  &YMW16ThinDisc::t2_b2);
