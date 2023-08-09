#include <cmath>
#include "hamunits.h"
#include "WMAP.h"

#include "helpers.h"


// ??????, implementation from Hammurabi (old)

std::array<double, 3>  WMAPMagneticField::at_position(const double &x, const double &y, const double &z) const { 
    
    std::array<double, 3>  B_vec3;
    double r = sqrt(x*x + y*y);

    if (r > b_r_max || r < b_r_min) { 
        return {0, 0, 0};
    }

	double phi = atan2(y, x);

    double psi_r = b_psi0*(M_PI/180.) + b_psi1*(M_PI/180.) * log(r/b_r0);
    double xsi_z = b_xsi0*(M_PI/180.) * tanh(z/b_z0);
    
    std::array<double, 3> B_cyl = {b_b0 * sin(psi_r) * cos(xsi_z), 
                                   b_b0 * cos(psi_r) * cos(xsi_z), 
                                   b_b0 * sin(xsi_z) };

    B_vec3 = Cyl2Cart(phi, B_cyl);
    


    // Antisymmetric, swap the signs.  The way my pitch angle is defined,
    // it seems this has to be swapped this way.  <------ hammurabi comment

    if (anti && z > 0) {
        B_vec3[0] *= (-1.); 
        B_vec3[1] *= (-1.); 
        B_vec3[2] *= (-1.);
        }
    return B_vec3;
}