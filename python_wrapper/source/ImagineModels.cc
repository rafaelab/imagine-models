#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>
#include <iostream>

//PYBIND11_MAKE_OPAQUE(std::vector<double>);
//PYBIND11_MAKE_OPAQUE(std::array<double, 3>);
//PYBIND11_MAKE_OPAQUE(std::array<int, 3>);

#include "include/regular_trampoline.h"
#include "../c_library/headers/Helix.h"
#include "../c_library/headers/Jaffe.h"
#include "../c_library/headers/Uniform.h"
#include "../c_library/headers/YMW.h"
#include "../c_library/headers/Sun.h"

if (HAVE_FFTW) {
  #include "include/random_trampoline.h"

  #include "../c_library/headers/EnsslinSteininger.h"
  #include "../c_library/headers/RandomJF12.h"
  #include "../c_library/headers/GaussianScalar.h"
  #include "../c_library/headers/LogNormal.h"
}


namespace py = pybind11;
using namespace pybind11::literals;

void FieldBases(py::module_ &);
void RegularFieldBases(py::module_ &);

void RegularJF12(py::module_ &);
void Helix(py::module_ &);
void Uniform(py::module_ &);
void Jaffe(py::module_ &);
void Sun2008(py::module_ &);


PYBIND11_MODULE(_ImagineModels, m) {
    m.doc() = "IMAGINE Model Library";

    FieldBases(m);
    RegularFieldBases(m);
    RegularJF12(m);
    Jaffe(m);
    Uniform(m);
    Sun2008(m);
    Helix(m);


/////////////////////////////////Regular Fields/////////////////////////////////


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


    py::class_<UniformMagneticField, RegularVectorField>(m, "UniformMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &UniformMagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("bx", &UniformMagneticField::bx)
        .def_readwrite("by", &UniformMagneticField::by)
        .def_readwrite("bz", &UniformMagneticField::bz);

    py::class_<Sun2008MagneticField, RegularVectorField>(m, "Sun2008MagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("at_position", &Sun2008MagneticField::at_position,  "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("b_B0", &Sun2008MagneticField::b_B0)
        .def_readwrite("b_Rsun", &Sun2008MagneticField::b_Rsun)
        .def_readwrite("b_R0", &Sun2008MagneticField::b_R0)
        .def_readwrite("b_z0", &Sun2008MagneticField::b_z0)
        .def_readwrite("b_Rc", &Sun2008MagneticField::b_Rc)
        .def_readwrite("b_Bc", &Sun2008MagneticField::b_Bc)
        .def_readwrite("b_pitch_deg", &Sun2008MagneticField::b_pitch_deg)

        .def_readwrite("bH_B0", &Sun2008MagneticField::bH_B0)
        .def_readwrite("bH_z0", &Sun2008MagneticField::bH_z0)
        .def_readwrite("bH_z1a", &Sun2008MagneticField::bH_z1a)
        .def_readwrite("bH_z1b", &Sun2008MagneticField::bH_z1b)
        .def_readwrite("bH_R0", &Sun2008MagneticField::bH_R0);


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

py::class_<YMW16, RegularScalarField>(m, "YMW16")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())

        .def("at_position", &YMW16::at_position, "x"_a, "y"_a, "z"_a, py::return_value_policy::move)

        .def_readwrite("r_warp", &YMW16::r_warp)
        .def_readwrite("r0", &YMW16::r0)
        .def_readwrite("t0_gamma_w", &YMW16::t0_gamma_w)
        // Thick disc
        .def_readwrite("t1_ad", &YMW16::t1_ad)
        .def_readwrite("t1_bd", &YMW16::t1_bd)
        .def_readwrite("t1_n1", &YMW16::t1_n1)
        .def_readwrite("t1_h1", &YMW16::t1_h1)
        // Thin disc
        .def_readwrite("t2_a2", &YMW16::t2_a2)
        .def_readwrite("t2_b2", &YMW16::t2_b2)
        .def_readwrite("t2_n2", &YMW16::t2_n2)
        .def_readwrite("t2_k2", &YMW16::t2_k2)
        // Spiral arms
        .def_readwrite("t3_b2s", &YMW16::t3_b2s)
        .def_readwrite("t3_ka", &YMW16::t3_ka)
        .def_readwrite("t3_aa", &YMW16::t3_aa)
        .def_readwrite("t3_ncn", &YMW16::t3_ncn)
        .def_readwrite("t3_wcn", &YMW16::t3_wcn)
        .def_readwrite("t3_thetacn", &YMW16::t3_thetacn)
        .def_readwrite("t3_nsg", &YMW16::t3_nsg)
        .def_readwrite("t3_wsg", &YMW16::t3_wsg)
        .def_readwrite("t3_thetasg", &YMW16::t3_thetasg)
        .def_readwrite("t3_rmin", &YMW16::t3_rmin)
        .def_readwrite("t3_phimin", &YMW16::t3_phimin)
        .def_readwrite("t3_tpitch", &YMW16::t3_tpitch)
        .def_readwrite("t3_cpitch", &YMW16::t3_cpitch)
        .def_readwrite("t3_narm", &YMW16::t3_narm)
        .def_readwrite("t3_warm", &YMW16::t3_warm)
        // Galactic center
        .def_readwrite("t4_ngc", &YMW16::t4_ngc)
        .def_readwrite("t4_agc", &YMW16::t4_agc)
        .def_readwrite("t4_hgc", &YMW16::t4_hgc)
        // Gum
        .def_readwrite("t5_kgn", &YMW16::t5_kgn)
        .def_readwrite("t5_ngn", &YMW16::t5_ngn)
        .def_readwrite("t5_wgn", &YMW16::t5_wgn)
        .def_readwrite("t5_agn", &YMW16::t5_agn)
        // Local bubble
        .def_readwrite("t6_j_lb", &YMW16::t6_j_lb)
        .def_readwrite("t6_nlb1", &YMW16::t6_nlb1)
        .def_readwrite("t6_detlb1", &YMW16::t6_detlb1)
        .def_readwrite("t6_wlb1", &YMW16::t6_wlb1)
        .def_readwrite("t6_hlb1", &YMW16::t6_hlb1)
        .def_readwrite("t6_thetalb1", &YMW16::t6_thetalb1)
        .def_readwrite("t6_nlb2", &YMW16::t6_nlb2)
        .def_readwrite("t6_detlb2", &YMW16::t6_detlb2)
        .def_readwrite("t6_wlb2", &YMW16::t6_wlb2)
        .def_readwrite("t6_hlb2", &YMW16::t6_hlb2)
        .def_readwrite("t6_thetalb2", &YMW16::t6_thetalb2)
        // Loop
        .def_readwrite("t7_nli", &YMW16::t7_nli)
        .def_readwrite("t7_rli", &YMW16::t7_rli)
        .def_readwrite("t7_wli", &YMW16::t7_wli)
        .def_readwrite("t7_detthetali", &YMW16::t7_detthetali)
        .def_readwrite("t7_thetali", &YMW16::t7_thetali);

if (HAVE_FFTW) {
/////////////////////////////////Random Fields/////////////////////////////////
// Random Vector Base Class
    py::class_<RandomVectorField, RandomField<std::array<double, 3>, std::array<double*, 3>>, PyRandomVectorField>(m, "RandomVectorField")
      .def(py::init<>())
      .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

      .def("on_grid", [](RandomVectorField &self, std::array<int, 3> &grid_shape,  std::array<double, 3>  &grid_reference_point, std::array<double, 3>  &grid_increment, int seed)  {
          std::array<double*, 3> f = self.on_grid(grid_shape, grid_reference_point, grid_increment, seed);
          size_t sx = grid_shape[0]; 
          size_t sy = grid_shape[1];
          size_t sz = grid_shape[2] + 1; // catches fftw zeropad (uneven)
          if (sz & 2) {
            sz += 1; // catches fftw zeropad (even)
          }
          auto lis = from_pointer_array_to_list_pyarray(f, sx, sy, sz, true);
          return lis;},
          py::kw_only(), py::arg("shape").noconvert(), py::arg("reference_point"), py::arg("increment"),  "seed"_a, 
          py::return_value_policy::take_ownership)

        .def("on_grid", [](RandomVectorField &self, int seed)  {
          std::array<double*, 3> f = self.on_grid(seed);
          size_t sx = self.shape[0];
          size_t sy = self.shape[1];
          size_t sz = self.shape[2] + 1; // catches fftw zeropad (uneven)
          if (sz & 2) {
            sz += 1; // catches fftw zeropad (even)
          }
          auto arr = from_pointer_array_to_list_pyarray(f, sx, sy, sz, true);
          return arr;}, 
          "seed"_a, 
          py::return_value_policy::take_ownership);

// Random Scalar Base Class
    py::class_<RandomScalarField, RandomField<double, double*>, PyRandomScalarField>(m, "RandomScalarField")
      .def(py::init<>())
      .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

      .def("on_grid", [](RandomScalarField &self, std::array<int, 3> &grid_shape,  std::array<double, 3>  &grid_reference_point, std::array<double, 3>  &grid_increment, int seed)  {
          double* f = self.on_grid(grid_shape, grid_reference_point, grid_increment, seed);
          size_t sx = grid_shape[0];
          size_t sy = grid_shape[1];
          size_t sz = grid_shape[2] + 1; // catches fftw zeropad (uneven)
          if (sz & 2) {
            sz += 1; // catches fftw zeropad (even)
          }
          std::cout<< "sz: " << sz << std::endl;
          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz, true);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;},
          py::kw_only(), py::arg("shape").noconvert(), py::arg("reference_point"), py::arg("increment"),  "seed"_a, 
          py::return_value_policy::take_ownership)

     .def("on_grid", [](RandomScalarField &self, int seed)  {
          double* f = self.on_grid(seed);
          size_t sx = self.shape[0];
          size_t sy = self.shape[1];
          size_t sz = self.shape[2] + 1; // catches fftw zeropad (uneven)
          if (sz & 2) {
            sz += 1; // catches fftw zeropad (even)
          }
          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz, true);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;}, 
          "seed"_a, 
          py::return_value_policy::take_ownership);
      

  /////////////////////////////////RandomFields/////////////////////////////////

    py::class_<JF12RandomField, RandomVectorField>(m, "JF12RandomField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readonly("regular_base", &JF12RandomField::regular_base)

        .def_readwrite("spectral_amplitude", &JF12RandomField::spectral_amplitude)
        .def_readwrite("spectral_offset", &JF12RandomField::spectral_offset)
        .def_readwrite("spectral_slope", &JF12RandomField::spectral_slope)

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
        .def_readwrite("Rmax", &JF12RandomField::Rmax)
        .def_readwrite("anisotropy_rho", &JF12RandomField::anisotropy_rho);


py::class_<ESRandomField, RandomVectorField>(m, "ESRandomField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("spectral_amplitude", &ESRandomField::spectral_amplitude)
        .def_readwrite("spectral_offset", &ESRandomField::spectral_offset)
        .def_readwrite("spectral_slope", &ESRandomField::spectral_slope)

        .def_readwrite("r0", &ESRandomField::r0)
        .def_readwrite("z0", &ESRandomField::z0);

/////////////////////////////Scalar Fields/////////////////////////////


py::class_<GaussianScalarField, RandomScalarField>(m, "GaussianScalarField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("mean", &GaussianScalarField::mean)
        .def_readwrite("spectral_amplitude", &GaussianScalarField::spectral_amplitude)
        .def_readwrite("spectral_offset", &GaussianScalarField::spectral_offset)
        .def_readwrite("spectral_slope", &GaussianScalarField::spectral_slope);


py::class_<LogNormalScalarField, RandomScalarField>(m, "LogNormalScalarField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("log_mean", &LogNormalScalarField::log_mean)
        .def_readwrite("spectral_amplitude", &LogNormalScalarField::spectral_amplitude)
        .def_readwrite("spectral_offset", &LogNormalScalarField::spectral_offset)
        .def_readwrite("spectral_slope", &LogNormalScalarField::spectral_slope);

  } // HAVE_FFTW
}