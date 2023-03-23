#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

class UniformMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        double bx = 0.;
        double by = 0.;
        double bz = 0.;

        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const {
            return std::array<double, 3> {bx, by, bz};
        }
 };
