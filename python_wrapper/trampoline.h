#include "../c_library/headers/hamunits.h"
#include "../c_library/headers/Field.h"
#include "../c_library/headers/RegularField.h"
#include "../c_library/headers/RandomField.h"
#include "../c_library/headers/RegularJF12.h"
#include "../c_library/headers/Helix.h"
#include "../c_library/headers/Jaffe.h"
//#include "../c_library/headers/ThermalElectronField.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <iostream>
// PYBIND11_MAKE_OPAQUE(std::vector<double>);

namespace py = pybind11;

// These classes are necessary to override virtual functions when binding abstract c++ classes
//FIELD CLASS MISSING

class PyVectorBaseField: public Field<std::array<double, 3>, std::array<double*, 3>> {
public:
    using Field<std::array<double, 3>, std::array<double*, 3>> ::Field;

    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::array<double, 3>, Field, at_position, x, y, z); }

    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) const override {PYBIND11_OVERRIDE_PURE(std::array<double*, 3>, Field, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) const override {PYBIND11_OVERRIDE_PURE(std::array<double*, 3>, Field, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

    std::array<double*, 3> on_grid(int seed) const override {PYBIND11_OVERRIDE_PURE(std::array<double*, 3>, Field, on_grid, seed); }
};


class PyScalarBaseField: public Field<double, double*> {
public:
    using Field<double, double*> ::Field;

    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::array<double, 3>, Field, at_position, x, y, z); }

    double* on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) const override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, grid_x, grid_y, grid_z, seed); }

    double* on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) const override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

    double* on_grid(int seed) const override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, seed); }
};


class PyRegularVectorField: public RegularVectorField {
public:
    using Field<std::array<double, 3>, std::array<double*, 3>>::Field; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::array<double, 3>, RegularVectorField, at_position, x, y, z); }

    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) const override {PYBIND11_OVERRIDE(std::array<double*, 3>, Field, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) const override {PYBIND11_OVERRIDE(std::array<double*, 3>, Field, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }
};

class PyRegularScalarField : public RegularScalarField {
public:
    using Field<py::array_t<double>, double>::Field; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, Field, at_position, x, y, z); }
    
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, Field, on_grid, grid_x, grid_y, grid_z); }
};

/*
class PyJF12MagneticField : public JF12MagneticField {
public:
    using JF12MagneticField::JF12MagneticField; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE(std::vector<double>, JF12MagneticField, at_position, x, y, z); }
    std::array<double*, 3> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE(std::vector<double>, JF12MagneticField, on_grid,  grid_x, grid_y, grid_z); }
};


class PyHelixMagneticField : public HelixMagneticField {
public:
    using HelixMagneticField::HelixMagneticField; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE(std::vector<double>, HelixMagneticField, at_position, x, y, z); }
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE(std::vector<double>, HelixMagneticField, on_grid,  grid_x, grid_y, grid_z); }
};

class PyJaffeMagneticField : public JaffeMagneticField {
public:
    using JaffeMagneticField::JaffeMagneticField; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE(std::vector<double>, JaffeMagneticField, at_position, x, y, z); }
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE(std::vector<double>, JaffeMagneticField, on_grid,  grid_x, grid_y, grid_z); }

};
*/
/*
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const {
      std::cout << "trampoline says hi " << std::endl;
      py::gil_scoped_acquire gil;
      py::function override = py::get_override(this, "on_grid");
      if (override) {
          py::array_t<double> obj = override(grid_x, grid_y, grid_z);
          std::cout << "trampline on grid obj size " << obj.size() << std::endl;
          return std::vector<double>(obj.data(), obj.data()+obj.size());
      }
      return on_grid(grid_x, grid_y, grid_z);
      }
*/
//class PyRandomField : public RandomField {
//public:
//    using RandomField::RandomField; // Inherit constructors
//    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, RandomField, at_position, x, y, z); }
//};