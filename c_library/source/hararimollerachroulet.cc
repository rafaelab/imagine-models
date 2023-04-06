#include <cmath>
#include "../headers/HarariMollerachRoulet.h"

#include "../headers/helpers.h"

std::array<double, 3>  HMRMagneticField::_at_position(const double &x, const double &y, const double &z, const HMRParams &p) const { 

   double r = sqrt(x*x + y*y);
   double phi = atan2(y, x);

   double f_z = (1. / (2.* cosh(z / p.b_z1))) + (1. / (2. * cosh(z / p.b_z2)));

   if (r < 0.0000000005) {
      r = 0.5;
      }

  double b_r = ( 3. * p.b_Rsun / r) * tanh(r / p.b_r1) * tanh(r / p.b_r1) * tanh(r / p.b_r1);

  double B_r_phi = b_r * cos(phi - ((1. / tan(p.b_p)) * log(r / p.b_epsilon0)));

  // B-field in cylindrical coordinates:
  std::array<double, 3> B_cyl{B_r_phi * sin(p.b_p) * f_z, 
                              B_r_phi * cos(p.b_p) * f_z , 
                              0.};

  std::array<double, 3> B_vec3{0, 0, 0};
  if (r <= p.b_r_max) {
    Cyl2Cart(phi, B_cyl, B_vec3);
  }
  
  return B_vec3;
}