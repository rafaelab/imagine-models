#include "../c_library/headers/MagneticField.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <iostream>
// PYBIND11_MAKE_OPAQUE(std::vector<double>);

namespace py = pybind11;

// These classes are necessary to override virtual functions when binding


class PyVectorField: public Field<py::array_t<double>, std::vector<double>> {
public:
    using Field<py::array_t<double>, std::vector<double>>::Field; // Inherit constructors
    std::vector<double> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, Field, at_position, x, y, z); }
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, Field, on_grid, grid_x, grid_y, grid_z); }
};

class PyScalarField : public Field<py::array_t<double>, double> {
public:
    using Field<py::array_t<double>, double>::Field; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, Field, at_position, x, y, z); }
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, Field, on_grid, grid_x, grid_y, grid_z); }
};

class PyRegularVectorField : public RegularField<py::array_t<double>, std::vector<double>> {
public:
    using RegularField<py::array_t<double>, std::vector<double>>::RegularField; // Inherit constructors
    std::vector<double> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, RegularField, at_position, x, y, z); }
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE(std::vector<double>, RegularField, on_grid,  grid_x, grid_y, grid_z); }
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
};

class PyRegularScalarField : public RegularField<py::array_t<double>, double> {
public:
    using RegularField<py::array_t<double>, double>::RegularField; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, RegularField, at_position, x, y, z); }
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE(std::vector<double>, RegularField, on_grid,  grid_x, grid_y, grid_z); }
};

//class PyRandomField : public RandomField {
//public:
//    using RandomField::RandomField; // Inherit constructors
//    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, RandomField, at_position, x, y, z); }
//};


class PyJF12MagneticField : public JF12MagneticField<py::array_t<double>> {
public:
    using JF12MagneticField<py::array_t<double>>::JF12MagneticField; // Inherit constructors
    std::vector<double> at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE(std::vector<double>, JF12MagneticField, at_position, x, y, z); }
    std::vector<double> on_grid(const py::array_t<double>& grid_x, const py::array_t<double>& grid_y, const py::array_t<double>& grid_z) const override {PYBIND11_OVERRIDE(std::vector<double>, JF12MagneticField, on_grid,  grid_x, grid_y, grid_z); }
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
