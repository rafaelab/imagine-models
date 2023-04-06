#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

struct UniformParams : Params {
    //parameter_names = {"ampx", "ampy"};

    double bx = 0.;
    double by = 0.;
    double bz = 0.;

    };



class UniformMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        UniformParams param;

        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const {
            return std::array<double, 3> {param.bx, param.by, param.bz};
        }
 };
