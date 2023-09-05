#include <cmath>
#include "hamunits.h"
#include "Helix.h"

vector HelixMagneticField::_at_position(const double &x, const double &y, const double &z, const HelixParams &p
) const {
  
  const double r{sqrt(x*x + y*y)}; // radius in cylindrical coordinates
  const double phi{atan2(y, x) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates#
 
  vector b {{0.0, 0.0, 0.0}};
  
  if ((r > p.rmin) && (r < p.rmax)) {
    b[0] = p.ampx * cos(phi); 
    b[1] = p.ampy * sin(phi); 
    b[2] = p.ampz;
    }
  return b;
};

