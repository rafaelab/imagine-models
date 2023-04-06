#include <cmath>
#include "../headers/HarariMollerachRoulet.h"

#include "../headers/helpers.h"

std::array<double, 3>  HMRMagneticField::at_position(const double &x, const double &y, const double &z) const { 

   double r = sqrt(x*x + y*y);
   double phi = atan2(y, x);

   double f_z = (1. / (2.* cosh(z / param.b_z1))) + (1. / (2. * cosh(z / param.b_z2)));

   if (r < 0.0000000005) {
      r = 0.5;
      }

  double b_r = ( 3. * param.b_Rsun / r) * tanh(r / param.b_r1) * tanh(r / param.b_r1) * tanh(r / param.b_r1);

  double B_r_phi = b_r * cos(phi - ((1. / tan(param.b_p)) * log(r / param.b_epsilon0)));

  // B-field in cylindrical coordinates:
  std::array<double, 3> B_cyl{B_r_phi * sin(param.b_p) * f_z, 
                              B_r_phi * cos(param.b_p) * f_z , 
                              0.};

  std::array<double, 3> B_vec3{0, 0, 0};
  if (r <= param.b_r_max) {
    Cyl2Cart(phi, B_cyl, B_vec3);
  }
  
  return B_vec3;
}