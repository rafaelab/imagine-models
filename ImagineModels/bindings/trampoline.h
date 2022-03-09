#include "../src/jf12.cc"

// These classes are necessary to override virtual functions

class PyMagneticField : public MagneticField {
public:
    using MagneticField::MagneticField; // Inherit constructors
};
class PyRegularMagneticField : public RegularMagneticField {
public:
    using RegularMagneticField::RegularMagneticField; // Inherit constructors
    std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(std::vector<double>, RegularMagneticField, evaluate_at_pos, pos); }
};

class PyJF12MagneticField : public JF12MagneticField {
public:
    using JF12MagneticField::JF12MagneticField; // Inherit constructors
    std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(std::vector<double>, JF12MagneticField, evaluate_at_pos, pos); }
};