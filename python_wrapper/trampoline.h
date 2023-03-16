#include "../c_library/headers/hamunits.h"
#include "../c_library/headers/Field.h"
#include "../c_library/headers/RegularField.h"
#include "../c_library/headers/RandomField.h"
#include "../c_library/headers/RegularJF12.h"
#include "../c_library/headers/EnsslinSteininger.h"
#include "../c_library/headers/RandomJF12.h"
#include "../c_library/headers/Helix.h"
#include "../c_library/headers/Jaffe.h"
//#include "../c_library/headers/ThermalElectronField.h"
#include <iostream>

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


    void allocate_memory(std::array<int, 3> shp, bool do_allocation, bool call_from_init) override {PYBIND11_OVERRIDE_PURE(void, Field, allocate_memory, shp, do_allocation, call_from_init); }
    
    void free_memory(bool do_deallocation) override {PYBIND11_OVERRIDE_PURE(void, Field, free_memory, do_deallocation); }
    
};


class PyScalarRandomFieldBase: public RandomField<double, double*> {
public:
    using RandomField<double, double*>:: RandomField; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, RandomField, at_position, x, y, z); }

    double* on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(double*, RandomField, on_grid, seed); }
    
    double* on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE_PURE(double*, RandomField, on_grid, grid_x, grid_y, grid_z, seed); }

    double* on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(double*, RandomField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

    double spatial_profile(const double &x, const double &y, const double &z) const override{PYBIND11_OVERRIDE_PURE(double, RandomField, spatial_profile, x, y, z); }

    void allocate_memory(std::array<int, 3> shp, bool do_allocation, bool call_from_init) override {PYBIND11_OVERRIDE_PURE(void, RandomField, allocate_memory, shp, do_allocation, call_from_init); }
    
    void free_memory(bool do_deallocation) override {PYBIND11_OVERRIDE_PURE(void, RandomField, free_memory, do_deallocation); }
    
};


class PyVectorFieldBase: public Field<std::array<double, 3>, std::array<double*, 3>> {
public:
    using Field<std::array<double, 3>, std::array<double*, 3>>:: Field; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(Array3Type, Field, at_position, x, y, z); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, seed); }
    
    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }


    void allocate_memory(std::array<int, 3> shp, bool do_allocation, bool call_from_init) override {PYBIND11_OVERRIDE_PURE(void, Field, allocate_memory, shp, do_allocation, call_from_init); }
    
    void free_memory(bool do_deallocation) override {PYBIND11_OVERRIDE_PURE(void, Field, free_memory, do_deallocation); }
    
};


class PyVectorRandomFieldBase: public RandomField<std::array<double, 3>, std::array<double*, 3>> {
public:
    using RandomField<std::array<double, 3>, std::array<double*, 3>>:: RandomField; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(Array3Type, Field, at_position, x, y, z); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, RandomField, on_grid, seed); }
    
    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, RandomField, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, RandomField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

    double spatial_profile(const double &x, const double &y, const double &z) const override{PYBIND11_OVERRIDE_PURE(double, RandomField, spatial_profile, x, y, z); }

    void allocate_memory(std::array<int, 3> shp, bool do_allocation, bool call_from_init) override {PYBIND11_OVERRIDE_PURE(void, RandomField, allocate_memory, shp, do_allocation, call_from_init); }
    
    void free_memory(bool do_deallocation) override {PYBIND11_OVERRIDE_PURE(void, RandomField, free_memory, do_deallocation); }
    
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


class PyRandomVectorField: public RandomVectorField {
public:
    using RandomVectorField:: RandomVectorField; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(Array3Type, RandomVectorField, at_position, x, y, z); }

    double spatial_profile(const double &x, const double &y, const double &z) const override{PYBIND11_OVERRIDE_PURE(double, RandomVectorField, spatial_profile, x, y, z); }

    void _on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(void, RandomVectorField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RandomVectorField, on_grid, seed); }
    
    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RandomVectorField, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RandomVectorField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }
};


class PyRandomScalarField: public RandomScalarField {
public:
    using RandomScalarField:: RandomScalarField; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, RandomScalarField, at_position, x, y, z); }

    double spatial_profile(const double &x, const double &y, const double &z) const override{PYBIND11_OVERRIDE_PURE(double, RandomScalarField, spatial_profile, x, y, z); }

    void _on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(void, RandomScalarField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

    double* on_grid(int seed) override {PYBIND11_OVERRIDE(double*, RandomScalarField, on_grid, seed); }
    
    double* on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE(double*, RandomScalarField, on_grid, grid_x, grid_y, grid_z, seed); }

    double* on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE(double*, RandomScalarField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

};

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
