#include "../c_library/headers/AbstractFields.h"

// These classes are necessary to override virtual functions when binding

class PyField : public Field {
public:
    using Field::Field; // Inherit constructors
};

class PyVectorField : public VectorField {
public:
    using VectorField::VectorField; // Inherit constructors
    std::vector<double> evaluate_model(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, VectorField, evaluate_model, x, y, z); }
};

class PyScalarField : public ScalarField {
public:
    using ScalarField::ScalarField; // Inherit constructors
    double evaluate_model(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, ScalarField, evaluate_model, x, y, z); }
};
