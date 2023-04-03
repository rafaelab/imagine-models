#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


// Fauvet magnetic field

class FauvetMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        double b_r_max = 20.; // kpc
        double b_r_min = 3.; // kpc
        double b_Rsun = 8.5; // kpc
        double b_b0 = 7.1; // muG
        double b_z0 = 1.; // kpc
        double b_r0 = 8.; // kpc
        double b_p = (M_PI/180.)*(-26.1); // degree
        double b_chi0 = (M_PI/180.)*22.4; // degree

        double h_b0 = 1.; // muG
        double h_z0 = 1.5; // kpc
        double h_z1a = .2; // kpc
        double h_z1b = .4; // kpc
        double h_r0 = 4.; // kpc

        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const;
 
};
