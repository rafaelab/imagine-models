#include <cmath>
#include "hamunits.h"
#include "HarariMollerachRoulet.h"

#include "helpers.h"

vector HMRMagneticField::_at_position(const double &x, const double &y, const double &z, const HMRParams &p) const { 

   double r = sqrt(x*x + y*y);
   double phi = atan2(y, x);

   double f_z = (1. / (2.* cosh(z / p.b_z1))) + (1. / (2. * cosh(z / p.b_z2)));

   if (r < 0.0000000005) {
      r = 0.5;
      }

  double b_r = ( 3. * p.b_Rsun / r) * tanh(r / p.b_r1) * tanh(r / p.b_r1) * tanh(r / p.b_r1);

  double B_r_phi = b_r * std::cos(phi - ((1. / std::tan(b_p*(M_PI/180.))) * std::log(r / b_epsilon0)));

  // B-field in cylindrical coordinates:
  std::array<double, 3> B_cyl{B_r_phi * std::sin(b_p*(M_PI/180.)) * f_z, 
                              B_r_phi * std::cos(b_p*(M_PI/180.)) * f_z , 
                              0.};

  std::array<double, 3> B_vec3;
  if (r <= b_r_max) {
    B_vec3 = Cyl2Cart(phi, B_cyl);
  }
  else {
    B_vec3 = {0, 0, 0};
  }
  
  return B_vec3;
}