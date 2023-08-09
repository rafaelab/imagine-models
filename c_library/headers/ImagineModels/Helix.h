#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


class HelixMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        double ampx = 0.;
        double ampy = 0.;
        double ampz = 0.;
        double rmax = 3.;
        double rmin = 0.;

        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const;
 };
