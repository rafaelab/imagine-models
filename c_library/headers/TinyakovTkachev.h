#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

struct TTParams : Params {
    //parameter_names = {"ampx", "ampy"};

    double b_b0 = 1.4; // muG
    double b_Rsun = 8.5; // kpc
    double b_r_max = 20.; // kpc
    double b_r_min = 4.; // kpc
    double b_d = -0.5; // kpc
    double b_z0 = 1.5; // kpc
    double b_p = (M_PI/180.)*(-8); // degree

    };



//Tinyakov and Tkachev (TT) https://arxiv.org/abs/astro-ph/0102101, implementation of https://arxiv.org/pdf/astro-ph/0510444.pdf (Kachelriess et al.)
class TTMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
        std::array<double, 3> _at_position(const double &x, const double &y, const double &z, const TTParams &p) const;
    public:
        using RegularVectorField :: RegularVectorField;

        TTParams param;

        std::array<double, 3> at_position(const double &x, const double &y, const double &z) const {
            return _at_position(x, y, z, this->param);
        }
 };
