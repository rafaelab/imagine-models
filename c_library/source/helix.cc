#include <cmath>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/MagneticField.h"

// Create helical magnetic field from rmin to rmax with form (bx cos(phi), by
// sin(phi), bz phi)
std::vector<double>  HelixMagneticField::evaluate_at_pos(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates
  const double phi{atan2(pos[1], pos[0]) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates#
  std::vector<double> b =  std::vector<double>{0.0, 0.0, 0.0};
  if ((r > rmin) && (r < rmax)) {
    b = std::vector<double>{bx * std::cos(phi), by * std::sin(phi), bz};
    }
  return b;
}

std::vector<double>  HelixMagneticField::_dbx_at_pos(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates
  const double phi{atan2(pos[1], pos[0]) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates

  std::vector<double> dbx =  std::vector<double>{0.0, 0.0, 0.0};
  if ((r > rmin) && (r < rmax)) {
    dbx = {std::cos(phi), 0., 0.};
  }
  return dbx;
}

std::vector<double>  HelixMagneticField::_dby_at_pos(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates
  const double phi{atan2(pos[1], pos[0]) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates

  std::vector<double> dby =  std::vector<double>{0.0, 0.0, 0.0};
  if ((r > rmin) && (r < rmax)) {
    dby = std::vector<double>{0., std::sin(phi), 0.};
  }
  return dby;
}

std::vector<double>  HelixMagneticField::_dbz_at_pos(const std::vector<double> &pos) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // radius in cylindrical coordinates

  std::vector<double> dbz =  std::vector<double>{0.0, 0.0, 0.0};
  if ((r > rmin) && (r < rmax)) {
    dbz = std::vector<double>{0., 0., 1.};
  }
  return dbz;
}
