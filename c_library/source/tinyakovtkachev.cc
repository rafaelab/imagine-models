#include <cmath>
#include "../headers/hamunits.h"
#include "../headers/TinyakovTkachev.h"

#include "../headers/helpers.h"

std::array<double, 3>  TTMagneticField::at_position(const double &x, const double &y, const double &z) const { 
    double r = std::sqrt(x*x + y*y);
    double phi = std::atan2(y, x);

    double beta = 1./std::tan(b_p*(M_PI/180.));
    if(std::abs(b_p*(M_PI/180.)) < 1.e-16) { 
        beta=1.;
        }

    double phase = (beta * log(1. + b_d / b_Rsun)) - M_PI/2.;
    //  double epsilon0=(b5_Rsun+b5_d)*exp(-(M_PI/2.)*tan(b5_p)); <-- hammurabi comment

    double sign;
    if (z < 0) { 
        sign=1.; 
    }
    else { 
        sign=-1.;
    }
    double f_z = sign * exp(-(abs(z) / b_z0));

    // there is a factor 1/cos(phase) difference between
    // the original TT and Kachelriess. <-- hammurabi comment
    //double b_r=b_b0*(b_Rsun/(r));
    double b_r = b_b0 * (b_Rsun / (r * cos(phase)));

    //if(r<b5_r_min) {b_r=b5_b0*(b5_Rsun/(b5_r_min));} <-- hammurabi comment
    if (r < b_r_min) {
        b_r = b_b0 * (b_Rsun / (b_r_min * cos(phase)));
        }
  
    double B_r_phi = b_r * cos(phi - beta*log(r/b_Rsun) + phase);
    //  if(r<1.e-26){B_r_phi = b_r*std::cos(phi+phase);} <-- hammurabi comment

    // B-field in cylindrical coordinates: <-- hammurabi comment
    std::array<double, 3> B_cyl{B_r_phi*std::sin(b_p*(M_PI/180.))*f_z, B_r_phi*std::cos(b_p*(M_PI/180.))*f_z , 0.};

    std::array<double, 3> B_vec3;
    if (r <= b_r_max) {
        B_vec3 = Cyl2Cart(phi, B_cyl);
    }
    else {
        B_vec3 = {0., 0., 0.};
    }

    return B_vec3;
    }
