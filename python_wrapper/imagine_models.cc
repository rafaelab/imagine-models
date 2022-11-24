//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
//#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <iostream>

#include "trampoline.h"
// #include "thermal_trampoline.h"


namespace py = pybind11;
using namespace pybind11::literals;


// helper function to avoid making a copy when returning a py::array_t
// author: https://github.com/YannickJadoul
// source: https://github.com/pybind/pybind11/issues/1042#issuecomment-642215028
template <typename Sequence>
inline py::array_t<typename Sequence::value_type> as_pyarray(Sequence &&seq) {
  auto size = seq.size();
  auto data = seq.data();
  std::unique_ptr<Sequence> seq_ptr =
      std::make_unique<Sequence>(std::move(seq));
  auto capsule = py::capsule(seq_ptr.get(), [](void *p) {
    std::unique_ptr<Sequence>(reinterpret_cast<Sequence *>(p));
  });
  seq_ptr.release();
  return py::array(size, data, capsule);
}


PYBIND11_MODULE(_ImagineModels, m) {
    m.doc() = "IMAGINE Model Library";

/////////////////////////////////Base classes/////////////////////////////////

// Vector Base Class
    py::class_<Field<py::array_t<double>, std::vector<double>>, PyVectorField>(m, "VectorField")
        .def(py::init<>());


// Scalar Base Class
    py::class_<Field<py::array_t<double>, double>, PyScalarField>(m, "ScalarField")
        .def(py::init<>());

/////////////////////////////////Regular Field Basc/////////////////////////////////

// Regular Vector Base Class
    py::class_<RegularField<py::array_t<double>, std::vector<double>>, Field<py::array_t<double>, std::vector<double>>, PyRegularVectorField>(m, "RegularVectorField")
        .def(py::init<>())
        // still produces a copy!
        .def("on_grid", [](RegularField<py::array_t<double>, std::vector<double>> &self, py::array_t<double> &grid_x, py::array_t<double> &grid_y, py::array_t<double> &grid_z)  {
          std::vector<double> f = self.on_grid(grid_x, grid_y, grid_z);
          std::cout << "on grid f size " << f.size() << std::endl;
          return as_pyarray(std::move(f));
        }, "grid_x"_a, "grid_y"_a, "grid_z"_a);

// Regular Scalar Base Class
    py::class_<RegularField<py::array_t<double>, double>, Field<py::array_t<double>, double>, PyRegularScalarField>(m, "RegularScalarField")
        .def(py::init<>())
        // still produces a copy!
        .def("on_grid", [](RegularField<py::array_t<double>, std::vector<double>> &self, py::array_t<double> &grid_x, py::array_t<double> &grid_y, py::array_t<double> &grid_z)  {
          std::vector<double> f = self.on_grid(grid_x, grid_y, grid_z);
          return as_pyarray(std::move(f));
        }, "grid_x"_a, "grid_y"_a, "grid_z"_a);
/////////////////////////////////Random Fields/////////////////////////////////

//TODO

/////////////////////////////////Regular Fields/////////////////////////////////

    py::class_<JF12MagneticField<py::array_t<double>>, RegularField<py::array_t<double>, std::vector<double>>, PyJF12MagneticField>(m, "JF12RegularField")
        .def(py::init<>())
        .def("at_position", &JF12MagneticField<py::array_t<double>>::at_position, "x"_a, "y"_a, "z"_a)
        .def_readwrite("b_arm_1", &JF12MagneticField<py::array_t<double>>::b_arm_1)
        .def_readwrite("b_arm_2", &JF12MagneticField<py::array_t<double>>::b_arm_2)
        .def_readwrite("b_arm_3", &JF12MagneticField<py::array_t<double>>::b_arm_3)
        .def_readwrite("b_arm_4", &JF12MagneticField<py::array_t<double>>::b_arm_4)
        .def_readwrite("b_arm_5", &JF12MagneticField<py::array_t<double>>::b_arm_5)
        .def_readwrite("b_arm_6", &JF12MagneticField<py::array_t<double>>::b_arm_6)
        .def_readwrite("b_arm_7", &JF12MagneticField<py::array_t<double>>::b_arm_7)
        .def_readwrite("b_ring", &JF12MagneticField<py::array_t<double>>::b_ring)
        .def_readwrite("h_disk", &JF12MagneticField<py::array_t<double>>::h_disk)
        .def_readwrite("w_disk", &JF12MagneticField<py::array_t<double>>::w_disk)
        .def_readwrite("Bn", &JF12MagneticField<py::array_t<double>>::Bn)
        .def_readwrite("Bs", &JF12MagneticField<py::array_t<double>>::Bs)
        .def_readwrite("rn", &JF12MagneticField<py::array_t<double>>::rn)
        .def_readwrite("rs", &JF12MagneticField<py::array_t<double>>::rs)
        .def_readwrite("wh", &JF12MagneticField<py::array_t<double>>::wh)
        .def_readwrite("z0", &JF12MagneticField<py::array_t<double>>::z0)
        .def_readwrite("B0_X", &JF12MagneticField<py::array_t<double>>::B0_X)
        .def_readwrite("Xtheta_const", &JF12MagneticField<py::array_t<double>>::Xtheta_const)
        .def_readwrite("rpc_X", &JF12MagneticField<py::array_t<double>>::rpc_X)
        .def_readwrite("r0_X", &JF12MagneticField<py::array_t<double>>::r0_X);
/*
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

    py::class_<RegularThermalElectronField, ThermalElectronField, PyRegularThermalElectronField>(m, "RegularThermalElectronField")
        .def(py::init<>())
        .def("_evaluate_grid", &RegularThermalElectronField::_evaluate_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a, "ev_at_pos"_a)
        .def("evaluate_grid", &RegularThermalElectronField::evaluate_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a);


    py::class_<YMW16Component, RegularThermalElectronField, PyYMW16Component>(m, "YMW16Component")
        .def(py::init<>())

        .def_readwrite("t1_ad",  &YMW16Component::t1_ad)
        .def_readwrite("t1_bd",  &YMW16Component::t1_bd);

    py::class_<YMW16ThickDisc, YMW16Component, PyYMW16ThickDisc>(m, "YMW16ThickDisc")
        .def(py::init<>())
        .def("evaluate_at_pos", &YMW16ThickDisc::evaluate_at_pos, "pos"_a)

        .def_readwrite("t1_n1",  &YMW16ThickDisc::t1_n1)
        .def_readwrite("t1_h1",  &YMW16ThickDisc::t1_h1);

   py::class_<YMW16ThinDisc, YMW16Component, PyYMW16ThinDisc>(m, "YMW16ThinDisc")
        .def(py::init<>())
        .def("evaluate_at_pos", &YMW16ThinDisc::evaluate_at_pos, "pos"_a)

        .def_readwrite("t2_n2",  &YMW16ThinDisc::t2_n2)
        .def_readwrite("t2_k2",  &YMW16ThinDisc::t2_k2)
        .def_readwrite("t2_a2",  &YMW16ThinDisc::t2_a2)
        .def_readwrite("t2_b2",  &YMW16ThinDisc::t2_b2);
*/
        }
