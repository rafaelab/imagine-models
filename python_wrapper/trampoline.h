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
using Array3Type = std::array<double, 3>;
using Array3PointerType = std::array<double*, 3>; // Only for PYBIND11_OVERRIDE_PURE macro, else gets confused by commas 

// These classes are necessary to override virtual functions when binding abstract c++ classes

class PyScalarFieldBase: public Field<double, double*> {
public:
    using Field<double, double*>:: Field; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, Field, at_position, x, y, z); }

    double* on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, seed); }
    
    double* on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, grid_x, grid_y, grid_z, seed); }

    double* on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }


    void allocate_memory(bool not_empty, int sz) override {PYBIND11_OVERRIDE_PURE(void, Field, allocate_memory, not_empty, sz); }
    
    void free_memory(bool not_empty) override {PYBIND11_OVERRIDE_PURE(void, Field, free_memory, not_empty); }
    
};


class PyVectorFieldBase: public Field<std::array<double, 3>, std::array<double*, 3>> {
public:
    using Field<std::array<double, 3>, std::array<double*, 3>>:: Field; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(Array3Type, Field, at_position, x, y, z); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, seed); }
    
    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }


    void allocate_memory(bool not_empty, int sz) override {PYBIND11_OVERRIDE_PURE(void, Field, allocate_memory, not_empty, sz); }
    
    void free_memory(bool not_empty) override {PYBIND11_OVERRIDE_PURE(void, Field, free_memory, not_empty); }
    
};


class PyRegularVectorField: public RegularVectorField {
public:
    using RegularVectorField:: RegularVectorField; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(Array3Type, RegularVectorField, at_position, x, y, z); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RegularVectorField, on_grid, seed); }
    
    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RegularVectorField, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RegularVectorField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }
    
};


class PyRegularScalarField : public RegularScalarField {
public:
    using RegularScalarField:: RegularScalarField; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, RegularScalarField, at_position, x, y, z); }
    
    double* on_grid(const int seed=0) override {PYBIND11_OVERRIDE_PURE(double*, RegularScalarField, on_grid,  seed); }

    double* on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, const int seed=0) override {PYBIND11_OVERRIDE_PURE(double*, RegularScalarField, on_grid, grid_x, grid_y, grid_z, seed); }

    double* on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed = 0) override {PYBIND11_OVERRIDE_PURE(double*, RegularScalarField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }
    
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