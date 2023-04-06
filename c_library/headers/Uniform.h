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
        vector  _at_position (const double &x, const double &y, const double &z, const UniformParams &p) const {
            return vector {{p.bx, p.by, p.bz}};
        }
    public:
        using RegularVectorField :: RegularVectorField;

        UniformParams param;

        vector at_position(const double &x, const double &y, const double &z) const {
            return _at_position(x, y, z, this->param);
        }
 };
