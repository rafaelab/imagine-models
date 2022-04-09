#include <cmath>
#include "../headers/hamunits.h"
#include "../headers/MagneticField.h"

// Create helical magnetic field from rmin to rmax with form (bx cos(phi), by
// sin(phi), bz phi)
std::vector<double>  HelixMagneticField::evaluate_at_pos(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates
  const double phi{atan2(pos[1], pos[0]) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates

  if ((r > rmin) && (r < rmax)) {
    return std::vector<double>{bx * std::cos(phi), by * std::sin(phi), bz};
  } else {
    return std::vector<double>{0.0, 0.0, 0.0};
  }
}

std::vector<double>  HelixMagneticField::dB_over_dbx(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates
  const double phi{atan2(pos[1], pos[0]) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates

  if ((r > rmin) && (r < rmax)) {
    return std::vector<double>{std::cos(phi), 0., 0.};
  } else {
    return std::vector<double>{0.0, 0.0, 0.0};
  }
}