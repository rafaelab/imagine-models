#include <ThermalElectronField.h>

// disk shape
// Parameters: t1_bd, t1_ad
double YMW16Component::gd(const double &rr) const {
  double gd{1.};
  if (rr > t1_bd) {
    gd = std::pow(
        1. / std::cosh((rr - t1_bd) / t1_ad),
        2);
  }
  return gd;
}

double YMW16Component::dgd_dbd(const double &rr) const {
  double gd{1.};
  if (rr > t1_bd) {
    double gd = Disc::gd(rr);

  }
  return 2*gd *std::tanh((rr - t1_bd) / t1_ad) / t1_ad ;
}

double YMW16Component::dgd_dad(const double &rr) const {
  double gd{1.};
  if (rr > t1_bd) {
    double gd = Disc::gd(rr);

  }
  return 2*gd * std::tanh((rr - t1_bd) / t1_ad) * (rr - t1_bd) / t1_ad / t1_ad;
}
