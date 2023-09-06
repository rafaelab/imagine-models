#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


// StanevBSS (HMR) see https://arxiv.org/abs/astro-ph/9607086

class StanevBSSMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        double b_p = -10; // degree
        double b_r_max = 20.; // kpc
        double b_r_min = 1.; // kpc
        double b_z0 = 1.; // kpc
        double b_r0 = 10.55; // kpc
        double b_Rsun = 8.5; // kpc
        double b_b0 = 6.; // muG
        double b_phi0 = M_PI; // radians

        vector at_position (const double &x, const double &y, const double &z) const;
 };