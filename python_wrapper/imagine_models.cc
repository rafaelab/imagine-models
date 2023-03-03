#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>
#include <iostream>

//PYBIND11_MAKE_OPAQUE(std::vector<double>);
//PYBIND11_MAKE_OPAQUE(std::array<double, 3>);
//PYBIND11_MAKE_OPAQUE(std::array<int, 3>);

#include "trampoline.h"
// #include "thermal_trampoline.h"


namespace py = pybind11;
using namespace pybind11::literals;


// helper functions to avoid making a copy when returning a py::array_t
// based on
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


inline py::list from_pointer_array_to_list_pyarray(std::array<double*, 3> seq, size_t arr_size_x, size_t arr_size_y, size_t arr_size_z) {
  size_t arr_size = arr_size_x*arr_size_y*arr_size_z;
  py::list li;
  for (int i = 0; i<3; ++i) {
    py::capsule free_when_done(seq[i], [](void *f) {
        std::unique_ptr<double>(reinterpret_cast<double*>(f));
        });
    py::array_t<double> arr = py::array(arr_size, seq[i], free_when_done);
    li.append(arr.reshape({arr_size_x, arr_size_y, arr_size_z}));
  }
  return li;
}

inline py::array_t<double> from_pointer_to_pyarray(double* data, size_t arr_size_x, size_t arr_size_y, size_t arr_size_z) {
    py::capsule free_when_done(data, [](void *f) {
        double *data = reinterpret_cast<double *>(f);
        std::cerr << "Element [0] = " << data[0] << "\n";
            std::cerr << "freeing memory @ " << f << "\n";
        delete[] data;
        });

  size_t arr_size = arr_size_x*arr_size_y*arr_size_z;
  py::array_t<double> arr = py::array(arr_size, data, free_when_done);
  return arr.reshape({arr_size_x, arr_size_y, arr_size_z});
}


PYBIND11_MODULE(_ImagineModels, m) {
    m.doc() = "IMAGINE Model Library";

////////////////////////////////Field Bases/////////////////////////////////

    py::class_<Field<std::array<double, 3>, std::array<double*, 3>>,  PyVectorFieldBase>(m, "VectorFieldBase");

    py::class_<Field<double, double*>,  PyScalarFieldBase>(m, "ScalarFieldBase");

    py::class_<RandomField<std::array<double, 3>, std::array<double*, 3>>,  PyVectorRandomFieldBase>(m, "VectorRandomFieldBase");

    py::class_<RandomField<double, double*>,  PyScalarRandomFieldBase>(m, "ScalarRandomFieldBase");

/////////////////////////////////Regular Field Bases/////////////////////////////////

// Regular Vector Base Class
    py::class_<RegularVectorField, Field<std::array<double, 3>, std::array<double*, 3>>, PyRegularVectorField>(m, "RegularVectorField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("on_grid", [](RegularVectorField &self, std::vector<double> &grid_x,  std::vector<double>  &grid_y, std::vector<double>  &grid_z)  {
            std::array<double*, 3> f = self.on_grid(grid_x, grid_y, grid_z);
            std::cout << "on grid f size " << f.size() << std::endl;
            int si = 3;
            size_t sx = grid_x.size();
            size_t sy = grid_y.size();
            size_t sz = grid_z.size();
            auto li = from_pointer_array_to_list_pyarray(f, sx, sy, sz);
            //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
            return li;},
            "grid_x"_a, "grid_y"_a, "grid_z"_a, py::return_value_policy::take_ownership)


        .def("on_grid", [](RegularVectorField &self, std::array<int, 3> &grid_shape,  std::array<double, 3>  &grid_zeropoint, std::array<double, 3>  &grid_increment)  {
          std::array<double*, 3> f = self.on_grid(grid_shape, grid_zeropoint, grid_increment);
          std::cout << "on grid f size " << f.size() << std::endl;
          int si = 3;
          size_t sx = grid_shape[0];
          size_t sy = grid_shape[1];
          size_t sz = grid_shape[2];
          auto arr = from_pointer_array_to_list_pyarray(std::move(f), sx, sy, sz);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;},
          "grid_shape"_a, "grid_zeropoint"_a, "grid_increment"_a, py::return_value_policy::take_ownership)

        .def("on_grid", [](RegularVectorField &self)  {
          std::array<double*, 3> f = self.on_grid();
          std::cout << "on grid f size " << f.size() << std::endl;
          int si = 3;
          size_t sx = self.shape[0];
          size_t sy = self.shape[1];
          size_t sz = self.shape[2];
          auto arr = from_pointer_array_to_list_pyarray(std::move(f), sx, sy, sz);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;}, py::return_value_policy::move);
        

// Regular Scalar Base Class
    py::class_<RegularScalarField, Field<double, double*>, PyRegularScalarField>(m, "RegularScalarField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("on_grid", [](RegularScalarField &self, std::vector<double> &grid_x,  std::vector<double>  &grid_y, std::vector<double>  &grid_z)  {
          double* f = self.on_grid(grid_x, grid_y, grid_z, 0);
          size_t sx = grid_x.size();
          size_t sy = grid_y.size();
          size_t sz = grid_z.size();
          auto arr = from_pointer_to_pyarray(f, sx, sy, sz);
          arr.resize({sx, sy, sz});
          return arr;},
          "grid_x"_a, "grid_y"_a, "grid_z"_a, py::return_value_policy::move)
        
        .def("on_grid", [](RegularScalarField &self, std::array<int, 3>  grid_shape, std::array<double, 3>  grid_zeropoint, std::array<double, 3>  grid_increment)  {
          double* f = self.on_grid(grid_shape, grid_zeropoint, grid_increment, 0);
          size_t sx = grid_shape[0];
          size_t sy = grid_shape[1];
          size_t sz = grid_shape[2];
          auto arr = from_pointer_to_pyarray(f, sx, sy, sz);
          return arr;},
          "grid_shape"_a, "grid_zeropoint"_a, "grid_increment"_a, py::return_value_policy::move)


        .def("on_grid", [](RegularScalarField &self)  {
          double* f = self.on_grid(0);
          size_t sx = self.shape[0];
          size_t sy = self.shape[1];
          size_t sz = self.shape[2];
          auto arr = from_pointer_to_pyarray(f, sx, sy ,sz);
          return arr;}, py::return_value_policy::move);

/////////////////////////////////Random Fields/////////////////////////////////

// Random Vector Base Class
    py::class_<RandomVectorField, RandomField<std::array<double, 3>, std::array<double*, 3>>, PyRandomVectorField>(m, "RandomVectorField")
      .def(py::init<>())
      .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

      .def("on_grid", [](RandomVectorField &self, std::array<int, 3> &grid_shape,  std::array<double, 3>  &grid_zeropoint, std::array<double, 3>  &grid_increment, int seed)  {
          std::array<double*, 3> f = self.on_grid(grid_shape, grid_zeropoint, grid_increment, seed);
          std::cout << "on grid f size " << f.size() << std::endl;
          int si = 3;
          size_t sx = grid_shape[0];
          size_t sy = grid_shape[1];
          size_t sz = grid_shape[2];
          auto arr = from_pointer_array_to_list_pyarray(std::move(f), sx, sy, sz);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;},
          "grid_shape"_a, "grid_zeropoint"_a, "grid_increment"_a, "seed"_a, py::return_value_policy::take_ownership)

        .def("on_grid", [](RandomVectorField &self, int seed)  {
          std::array<double*, 3> f = self.on_grid(seed);
          std::cout << "on grid f size " << f.size() << std::endl;
          int si = 3;
          size_t sx = self.shape[0];
          size_t sy = self.shape[1];
          size_t sz = self.shape[2];
          auto arr = from_pointer_array_to_list_pyarray(std::move(f), sx, sy, sz);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;}, "seed"_a, py::return_value_policy::move);

// Random Scalar Base Class
    py::class_<RandomScalarField, RandomField<double, double*>, PyRandomScalarField>(m, "RandomScalarField")
      .def(py::init<>())
      .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

      .def("on_grid", [](RandomScalarField &self, std::array<int, 3> &grid_shape,  std::array<double, 3>  &grid_zeropoint, std::array<double, 3>  &grid_increment, int seed)  {
          double* f = self.on_grid(grid_shape, grid_zeropoint, grid_increment, seed);
          size_t sx = grid_shape[0];
          size_t sy = grid_shape[1];
          size_t sz = grid_shape[2];
          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;},
          "grid_shape"_a, "grid_zeropoint"_a, "grid_increment"_a, "seed"_a, py::return_value_policy::take_ownership)

     .def("on_grid", [](RandomScalarField &self, int seed)  {
          double* f = self.on_grid(seed);
          size_t sx = self.shape[0];
          size_t sy = self.shape[1];
          size_t sz = self.shape[2];
          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;}, "seed"_a, py::return_value_policy::move);
      

/////////////////////////////////Regular Fields/////////////////////////////////

    py::class_<JF12MagneticField, RegularVectorField>(m, "JF12RegularField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &JF12MagneticField::at_position, "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
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

    py::class_<HelixMagneticField, RegularVectorField>(m, "HelixMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &HelixMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
  //      .def("dampx_grid", &HelixMagneticField::dampx_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a)
  //      .def("dampy_grid", &HelixMagneticField::dampy_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a)
  //      .def("dampz_grid", &HelixMagneticField::dampz_grid, "grid_x"_a, "grid_y"_a, "grid_z"_a)
        .def_readwrite("ampx", &HelixMagneticField::ampx)
        .def_readwrite("ampy", &HelixMagneticField::ampy)
        .def_readwrite("ampz", &HelixMagneticField::ampz)
        .def_readwrite("rmin", &HelixMagneticField::rmin)
        .def_readwrite("rmax", &HelixMagneticField::rmax);

    py::class_<JaffeMagneticField, RegularVectorField>(m, "JaffeMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &JaffeMagneticField::at_position, "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
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



  /////////////////////////////////RandomFields/////////////////////////////////

    py::class_<JF12RandomField, RandomVectorField>(m, "JF12RandomField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("rms", &JF12RandomField::rms)
        .def_readwrite("k0", &JF12RandomField::k0)
        .def_readwrite("k1", &JF12RandomField::k1)
        .def_readwrite("a0", &JF12RandomField::a0)
        .def_readwrite("a1", &JF12RandomField::a1)
        .def_readwrite("rho", &JF12RandomField::rho)

        .def_readwrite("b0_1", &JF12RandomField::b0_1)
        .def_readwrite("b0_2", &JF12RandomField::b0_2)
        .def_readwrite("b0_3", &JF12RandomField::b0_3)
        .def_readwrite("b0_4", &JF12RandomField::b0_4)
        .def_readwrite("b0_5", &JF12RandomField::b0_5)
        .def_readwrite("b0_6", &JF12RandomField::b0_6)
        .def_readwrite("b0_7", &JF12RandomField::b0_7)
        .def_readwrite("b0_8", &JF12RandomField::b0_8)
        .def_readwrite("b0_int", &JF12RandomField::b0_int)
        .def_readwrite("b0_halo", &JF12RandomField::b0_halo)
        .def_readwrite("r0_halo", &JF12RandomField::r0_halo)
        .def_readwrite("z0_halo", &JF12RandomField::z0_halo)
        .def_readwrite("z0_spiral", &JF12RandomField::z0_spiral)
        .def_readwrite("rho_GC", &JF12RandomField::rho_GC)
        .def_readwrite("Rmax", &JF12RandomField::Rmax);


py::class_<ESRandomField, RandomVectorField>(m, "ESRandomField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("rms", &ESRandomField::rms)
        .def_readwrite("k0", &ESRandomField::k0)
        .def_readwrite("k1", &ESRandomField::k1)
        .def_readwrite("a0", &ESRandomField::a0)
        .def_readwrite("a1", &ESRandomField::a1)
        .def_readwrite("rho", &ESRandomField::rho)

        .def_readwrite("r0", &ESRandomField::r0)
        .def_readwrite("z0", &ESRandomField::z0);
/////////////////////////////Thermal Electron Field/////////////////////////////
/*
    py::class_<YMW16, RegularField<py::array_t<double>, double>, PyYMW16>(m, "YMW16Component")
        .def(py::init<>())
        .def("at_position", &YMW16ThickDisc::evaluate_at_pos, "pos"_a)

        .def_readwrite("t1_ad",  &YMW16Component::t1_ad)
        .def_readwrite("t1_bd",  &YMW16Component::t1_bd);

    py::class_<YMW16ThickDisc, YMW16Component, PyYMW16ThickDisc>(m, "YMW16ThickDisc")
        .def(py::init<>())


        .def_readwrite("t1_n1",  &YMW16ThickDisc::t1_n1)
        .def_readwrite("t1_h1",  &YMW16ThickDisc::t1_h1);

   py::class_<YMW16ThinDisc, YMW16Component, PyYMW16ThinDisc>(m, "YMW16ThinDisc")
        .def(py::init<>())
        .def("at_posiiton", &YMW16ThinDisc::evaluate_at_pos, "pos"_a)

        .def_readwrite("t2_n2",  &YMW16ThinDisc::t2_n2)
        .def_readwrite("t2_k2",  &YMW16ThinDisc::t2_k2)
        .def_readwrite("t2_a2",  &YMW16ThinDisc::t2_a2)
        .def_readwrite("t2_b2",  &YMW16ThinDisc::t2_b2);
*/
        }
