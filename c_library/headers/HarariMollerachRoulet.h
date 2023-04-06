#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


// Harari, Mollerach, Roulet (HMR) see https://arxiv.org/abs/astro-ph/9906309, implementation of https://arxiv.org/pdf/astro-ph/0510444.pdf 

struct HMRParams : Params {
    //parameter_names = {"ampx", "ampy"};

    double b_Rsun = 8.5; // kpc
    double b_r_max = 20.; // kpc
    double b_z1 = 0.3; // kpc
    double b_z2 = 4.; // kpc
    double b_r1 = 2.; // kpc
    double b_p = (M_PI/180.)*(-10); // degree
    double b_epsilon0 = 10.55; // kpc
    };



class HMRMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        HMRParams param;

        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const;
 };