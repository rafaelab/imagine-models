#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


// WMAP magnetic field

class WMAPMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        double b_r_max = 20.; // kpc
        double b_r_min = 3.; // kpc
        double b_Rsun = 8.5; // kpc
        double b_b0 = 6.; // kpc
        double b_z0 = 1.; // kpc
        double b_r0 = 8.; // kpc
        double b_psi0 = (M_PI/180.)*27; // degree
        double b_psi1 = (M_PI/180.)*0.9; // degree
        double b_xsi0 = (M_PI/180.)*25; // degree

        bool anti = false;


        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const;
 };
