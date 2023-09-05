#include <cmath>
#include "hamunits.h"
#include "Fauvet.h"

#include "helpers.h"


// ??????, implementation from Hammurabi (old)

std::array<double, 3>  FauvetMagneticField::at_position(const double &x, const double &y, const double &z) const { 
    
    std::array<double, 3>  B_vec3;
    double r = sqrt(x*x + y*y);

    if (r > b_r_max || r < b_r_min) { 
        return {0, 0, 0};
    }

	double phi = atan2(y, x);
    double chi_z = b_chi0*(M_PI/180.) * tanh(z/b_z0);
    double beta=1./tan(b_p*(M_PI/180.));

    // B-field in cylindrical coordinates:
    std::array<double, 3> B_cyl{b_b0 * cos(phi + beta * log(r/b_r0)) * sin(b_p*(M_PI/180.)) * cos(chi_z),
                               -b_b0 * cos(phi + beta * log(r/b_r0)) * cos(b_p*(M_PI/180.)) * cos(chi_z), 
                                b_b0 * sin(chi_z)};
  

    // Taking into account the halo field
    double h_z1;
    if (abs(z) < h_z0) {
        h_z1 = h_z1a;
    }
    else { 
        h_z1 = h_z1b;
    }
        
    double hf_piece1 = (h_z1 * h_z1) / (h_z1 * h_z1 + (abs(z) - h_z0) * (abs(z) - h_z0));
    double hf_piece2 = exp(-(r-h_r0)/(h_r0));

    double halo_field = h_b0 * hf_piece1 * (r / b_r0) * hf_piece2;
    B_cyl[1] += halo_field;

    B_vec3 = Cyl2Cart(phi, B_cyl);
    return B_vec3;
};

