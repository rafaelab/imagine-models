#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


struct TTParams : Params {

    number b_b0 = 1.4; // muG
    number b_Rsun = 8.5; // kpc
    double b_r_max = 20.; // kpc
    double b_r_min = 4.; // kpc
    number b_d = -0.5; // kpc
    number b_z0 = 1.5; // kpc
    number b_p = -8; // degree

    };


//Tinyakov and Tkachev (TT) https://arxiv.org/abs/astro-ph/0102101, implementation of https://arxiv.org/pdf/astro-ph/0510444.pdf (Kachelriess et al.)
class TTMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
        vector  _at_position(const double &x, const double &y, const double &z, const TTParams &p) const;
    public:
        using RegularVectorField :: RegularVectorField;

        TTParams param;

        vector at_position(const double &x, const double &y, const double &z) const {
            return _at_position(x, y, z, this->param);
        }
        

        #if autodiff_FOUND
            Eigen::MatrixXd _derivative(const double &x, const double &y, const double &z,  TTParams &p) {
                vector out;
                Eigen::MatrixXd deriv = ad::jacobian([&](auto _x, auto _y, auto _z, auto _p) {return this->_at_position(_x, _y, _z, _p);}, ad::wrt(p.b_b0, p.b_d, p.b_p, p.b_Rsun, p.b_z0), ad::at(x, y, z, p), out);  
                return deriv;
            }

            Eigen::MatrixXd derivative(const double &x, const double &y, const double &z) {
                return _derivative(x, y, z, this->param);
            }
        #endif       
 };
