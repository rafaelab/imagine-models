#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


//Sun et al. A&A V.477 2008 ASS+RING model magnetic field
class Sun2008MagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;
        double b_B0;
        // Note, redundant with SunPosition. However
        // might want to tune it independently anyways.
        double b_Rsun; 
        double b_R0;
        double b_z0;
        double b_Rc;
        double b_Bc;
        double b_pitch_deg;
        
        double bH_B0;
        double bH_z0;
        double bH_z1a;
        double bH_z1b;
        double bH_R0;

        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const;
 };
