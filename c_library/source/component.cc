#include <cmath>

#include "../headers/ThermalElectronField.h"

// disk shape
// Parameters: t1_bd, t1_ad
double YMW16Component::gd(const double &rr) const {
  double gdv{1.};
  if (rr > t1_bd) {
    gdv = std::pow(
        1. / std::cosh((rr - t1_bd) / t1_ad),
        2);
  }
  return gdv;
}

double YMW16Component::dgd_dbd(const double &rr) const {
  double gdv{1.};
  if (rr > t1_bd) {
    double gdv = gd(rr);

  }
  return 2 * gdv *std::tanh((rr - t1_bd) / t1_ad) / t1_ad ;
}

double YMW16Component::dgd_dad(const double &rr) const {
  double gdv{1.};
  if (rr > t1_bd) {
    double gdv = gd(rr);

  }
  return 2 * gdv * std::tanh((rr - t1_bd) / t1_ad) * (rr - t1_bd) / t1_ad / t1_ad;
}
