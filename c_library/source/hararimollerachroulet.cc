#include <cmath>
#include "hamunits.h"
#include "HarariMollerachRoulet.h"

#include "helpers.h"

vector HMRMagneticField::_at_position(const double &x, const double &y, const double &z, const HMRParams &p) const { 

  vector B_vec3{{0, 0, 0}};

  double r = std::sqrt(x*x + y*y);
  double phi = std::atan2(y, x);

  number f_z = (1. / (2.* cosh(z / p.b_z1))) + (1. / (2. * cosh(z / p.b_z2)));

  if (r < 0.0000000005) {
    r = 0.5;
    }

  number b_r = ( 3. * p.b_Rsun / r) * tanh(r / p.b_r1) * tanh(r / p.b_r1) * tanh(r / p.b_r1);

  number B_r_phi = b_r * cos(phi - ((1. / tan(p.b_p*(M_PI/180.))) * log(r / p.b_epsilon0)));

  // B-field in cylindrical coordinates:
  vector B_cyl{{B_r_phi * sin(p.b_p*(M_PI/180.)) * f_z, 
               B_r_phi * cos(p.b_p*(M_PI/180.)) * f_z , 
               0.}
               };

  B_vec3 = Cyl2Cart(phi, B_cyl);
  
  return B_vec3;
}