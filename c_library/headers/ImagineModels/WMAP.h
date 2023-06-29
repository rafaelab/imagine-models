#include <functional>
#include <cmath>

#include "ImagineModels/Field.h"
#include "ImagineModels/RegularField.h"


// WMAP magnetic field

class WMAPMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        double b_r_max = 20.; // kpc
        double b_r_min = 3.; // kpc
        double b_Rsun = 8.5; // kpc
        double b_b0 = 6.; // muG
        double b_z0 = 1.; // kpc
        double b_r0 = 8.; // kpc
        double b_psi0 = 27; // degree
        double b_psi1 = 0.9; // degree
        double b_xsi0 = 25; // degree

        bool anti = false;


        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const;
 };
