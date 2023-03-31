#include "../c_library/headers/hamunits.h"
#include "../c_library/headers/Field.h"

#include "../c_library/headers/RandomField.h"


#include <iostream>

using Array3Type = std::array<double, 3>;
using Array3PointerType = std::array<double*, 3>; // Only for PYBIND11_OVERRIDE_PURE macro, else gets confused by commas 

// These classes are necessary to override virtual functions when binding abstract c++ classes


class PyScalarRandomFieldBase: public RandomField<double, double*> {
public:
    using RandomField<double, double*>:: RandomField; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, RandomField, at_position, x, y, z); }

    double* on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(double*, RandomField, on_grid, seed); }
    
    double* on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE(double*, RandomField, on_grid, grid_x, grid_y, grid_z, seed); }

    double* on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(double*, RandomField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

    double spatial_profile(const double &x, const double &y, const double &z) const override{PYBIND11_OVERRIDE_PURE(double, RandomField, spatial_profile, x, y, z); }

    double calculate_fourier_sigma(const double &abs_k) const override{PYBIND11_OVERRIDE_PURE(double, RandomField, calculate_fourier_sigma, abs_k); }

    double* allocate_memory(std::array<int, 3> shp) override {PYBIND11_OVERRIDE_PURE(double*, RandomField, allocate_memory, shp); }
    
    void free_memory(double* grid_eval) override {PYBIND11_OVERRIDE_PURE(void, RandomField, free_memory, grid_eval); }
    
};


class PyVectorRandomFieldBase: public RandomField<std::array<double, 3>, std::array<double*, 3>> {
public:
    using RandomField<std::array<double, 3>, std::array<double*, 3>>:: RandomField; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(Array3Type, Field, at_position, x, y, z); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, RandomField, on_grid, seed); }
    
    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, RandomField, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, RandomField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

    double spatial_profile(const double &x, const double &y, const double &z) const override{PYBIND11_OVERRIDE_PURE(double, RandomField, spatial_profile, x, y, z); }

    double calculate_fourier_sigma(const double &abs_k) const override{PYBIND11_OVERRIDE_PURE(double, RandomField, calculate_fourier_sigma, abs_k); }

    std::array<double*, 3> allocate_memory(std::array<int, 3> shp) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, RandomField, allocate_memory, shp); }
    
    void free_memory(std::array<double*, 3> grid_eval) override {PYBIND11_OVERRIDE_PURE(void, RandomField, free_memory, grid_eval); }
    
};



class PyRandomVectorField: public RandomVectorField {
public:
    using RandomVectorField:: RandomVectorField; // Inherit constructors
    std::array<double, 3> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(Array3Type, RandomVectorField, at_position, x, y, z); }

    double spatial_profile(const double &x, const double &y, const double &z) const override{PYBIND11_OVERRIDE_PURE(double, RandomVectorField, spatial_profile, x, y, z); }

    double calculate_fourier_sigma(const double &abs_k) const override{PYBIND11_OVERRIDE_PURE(double, RandomVectorField, calculate_fourier_sigma, abs_k); }

    void _on_grid(std::array<double*, 3> grid_eval, const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, const int seed) override {PYBIND11_OVERRIDE_PURE(void, RandomVectorField, on_grid, grid_eval, grid_shape, grid_zeropoint, grid_increment, seed); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RandomVectorField, on_grid, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, const int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RandomVectorField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }
};


class PyRandomScalarField: public RandomScalarField {
public:
    using RandomScalarField:: RandomScalarField; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, RandomScalarField, at_position, x, y, z); }

    double spatial_profile(const double &x, const double &y, const double &z) const override{PYBIND11_OVERRIDE_PURE(double, RandomScalarField, spatial_profile, x, y, z); }

    double calculate_fourier_sigma(const double &abs_k) const override{PYBIND11_OVERRIDE_PURE(double, RandomScalarField, calculate_fourier_sigma, abs_k); }

    void _on_grid(double* grid_eval, const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, const int seed) override {PYBIND11_OVERRIDE_PURE(void, RandomScalarField, on_grid,  grid_eval, grid_shape, grid_zeropoint, grid_increment, seed); }

    double* on_grid(const int seed) override {PYBIND11_OVERRIDE(double*, RandomScalarField, on_grid, seed); }

    double* on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, const int seed) override {PYBIND11_OVERRIDE(double*, RandomScalarField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }

};

