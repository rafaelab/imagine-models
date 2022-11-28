#include <cmath>
#include "../headers/hamunits.h"
#include "../headers/MagneticField.h"

// Create helical magnetic field from rmin to rmax with form (bx cos(phi), by
// sin(phi), bz phi)
double*  HelixMagneticField::evaluate_model(const double &x, const double &y, const double &z) const {
  const double r{sqrt(x*x] + y*y)}; // radius in cylindrical coordinates
  const double phi{atan2(y, x) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates#
  std::vector<double> b =  std::vector<double>{0.0, 0.0, 0.0};
  if ((r > rmin) && (r < rmax)) {
    b = std::vector<double>{ampx * std::cos(phi), ampy * std::sin(phi), ampz};
    }
  return b;
};
/*
std::vector<double>  HelixMagneticField::_dampx_at_pos(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates
  const double phi{atan2(pos[1], pos[0]) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates

  std::vector<double> dampx =  std::vector<double>{0.0, 0.0, 0.0};
  if ((r > rmin) && (r < rmax)) {
    dampx = std::vector<double>{std::cos(phi), 0., 0.};
  }
  return dampx;
}

std::vector<double>  HelixMagneticField::_dampy_at_pos(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates
  const double phi{atan2(pos[1], pos[0]) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates

  std::vector<double> dampy =  std::vector<double>{0.0, 0.0, 0.0};
  if ((r > rmin) && (r < rmax)) {
    dampy = std::vector<double>{0., std::sin(phi), 0.};
  }
  return dampy;
}

std::vector<double>  HelixMagneticField::_dampz_at_pos(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates

  std::vector<double> dampz =  std::vector<double>{0.0, 0.0, 0.0};
  if ((r > rmin) && (r < rmax)) {
    dampz = std::vector<double>{0., 0., 1.};
  }
  return dampz;
}
*/
