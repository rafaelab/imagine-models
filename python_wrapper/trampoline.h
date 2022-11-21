#include "../c_library/headers/MagneticField.h"

// These classes are necessary to override virtual functions when binding

class PyVectorField : public Field<std::vector<double>, std::vector<double>> {
public:
    using Field<std::vector<double>, std::vector<double>>::Field; // Inherit constructors
    std::vector<double> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, Field, at_position, x, y, z); }
    std::vector<double> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, Field, on_grid, grid_x, grid_y, grid_z); }
};

class PyScalarField : public Field<std::vector<double>, double> {
public:
    using Field<std::vector<double>, double>::Field; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, Field, at_position, x, y, z); }
    std::vector<double> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, Field, on_grid, grid_x, grid_y, grid_z); }
};


class PyRegularVectorField : public RegularField<std::vector<double>, std::vector<double>> {
public:
    using RegularField<std::vector<double>, std::vector<double>>::RegularField; // Inherit constructors
    std::vector<double> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, RegularField, at_position, x, y, z); }
    std::vector<double> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z) const override {PYBIND11_OVERRIDE(std::vector<double>, RegularField, on_grid,  grid_x, grid_y, grid_z); }
};

class PyRegularScalarField : public RegularField<std::vector<double>, double> {
public:
    using RegularField<std::vector<double>, double>::RegularField; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, RegularField, at_position, x, y, z); }
    std::vector<double> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z) const override {PYBIND11_OVERRIDE(std::vector<double>, RegularField, on_grid,  grid_x, grid_y, grid_z); }
};

//class PyRandomField : public RandomField {
//public:
//    using RandomField::RandomField; // Inherit constructors
//    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, RandomField, at_position, x, y, z); }
//};


class PyJF12MagneticField : public JF12MagneticField<std::vector<double>> {
public:
    using JF12MagneticField<std::vector<double>>::JF12MagneticField; // Inherit constructors
    std::vector<double> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE(std::vector<double>, JF12MagneticField, at_position, x, y, z); }
};

/*
class PyHelixMagneticField : public HelixMagneticField {
public:
    using HelixMagneticField::HelixMagneticField; // Inherit constructors
    std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(std::vector<double>, HelixMagneticField, evaluate_at_pos, pos); }

};

class PyJaffeMagneticField : public JaffeMagneticField {
public:
    using JaffeMagneticField::JaffeMagneticField; // Inherit constructors
    std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(std::vector<double>, JaffeMagneticField, evaluate_at_pos, pos); }

};
*/
