#include <cmath>
#include "../headers/hamunits.h"
#include "../headers/StanevBSS.h"

#include "../headers/helpers.h"


// https://arxiv.org/abs/astro-ph/9607086, implementation from Hammurabi (old)
std::array<double, 3>  StanevBSSMagneticField::at_position(const double &x, const double &y, const double &z) const { 
    
    std::array<double, 3>  B_vec3;
    double r = sqrt(x*x + y*y);
	double phi = atan2(y, x);

    if (r > b_r_max || r < b_r_min) { 
        return {0, 0, 0};
    }

    double phi_prime = b_phi0 - phi;  // PHIprime running clock-wise from neg. x-axis
    double beta = 1. / tan(b_p);

    double B_0 = b_b0 * b_Rsun / 4. ; 
    if (r > 4.) {
        B_0 = b_b0 * b_Rsun /r;}

    std::array<double, 3> B_cyl = {B_0 * cos(phi_prime -beta * log(r/b_r0)) * sin(b_p) * exp(-std::abs(z)/b_z0),
				                  -B_0 * cos(phi_prime -beta * log(r/b_r0)) * cos(b_p) * exp(-std::abs(z)/b_z0),
				                   0.};

    B_vec3 = Cyl2Cart(phi, B_cyl);
    return B_vec3;
}