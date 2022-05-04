#include "../c_library/headers/ThermalElectronField.h"

// These classes are necessary to override virtual functions when binding

class PyThermalElectronField : public ThermalElectronField{
public:
    using ThermalElectronField::ThermalElectronField; // Inherit constructors
};
class PyRegularThermalElectronField : public RegularThermalElectronField {
public:
    using RegularThermalElectronField::RegularThermalElectronField; // Inherit constructors
    double evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE_PURE(double, ThermalElectronField, evaluate_at_pos, pos); }
};

class PyYMW16ThinDisc : public YMW16ThinDisc {
public:
    using YMW16ThinDisc::YMW16ThinDisc; // Inherit constructors
    double evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(double, YMW16ThinDisc, evaluate_at_pos, pos); }
};

class PyYMW16ThickDisc : public YMW16ThickDisc {
public:
    using YMW16ThickDisc::YMW16ThickDisc; // Inherit constructors
    double evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(double, YMW16ThickDisc, evaluate_at_pos, pos); }

};
