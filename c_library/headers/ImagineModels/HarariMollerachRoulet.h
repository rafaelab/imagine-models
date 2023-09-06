#ifndef HMR_H
#define HMR_H

#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


// Harari, Mollerach, Roulet (HMR) see https://arxiv.org/abs/astro-ph/9906309, implementation of https://arxiv.org/pdf/astro-ph/0510444.pdf 

struct HMRParams : Params {

    number b_Rsun = 8.5; // kpc
    double b_r_max = 20.; // kpc
    number b_z1 = 0.3; // kpc
    number b_z2 = 4.; // kpc
    number b_r1 = 2.; // kpc
    number b_p = -10; // degree
    number b_epsilon0 = 10.55; // kpc

    };

class HMRMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
        vector  _at_position(const double &x, const double &y, const double &z, const HMRParams &p) const;
    public:
        using RegularVectorField :: RegularVectorField;

        HMRParams param;

        vector at_position(const double &x, const double &y, const double &z) const {
            return _at_position(x, y, z, this->param);
        }

        #if autodiff_FOUND
            Eigen::MatrixXd _derivative(const double &x, const double &y, const double &z,  HMRParams &p) {
                vector out;
                Eigen::MatrixXd deriv = ad::jacobian([&](auto _x, auto _y, auto _z, auto _p) {return this->_at_position(_x, _y, _z, _p);}, ad::wrt(p.b_Rsun, p.b_p, p.b_r1, p.b_z1, p.b_z2, p.b_epsilon0), ad::at(x, y, z, p), out);  
                return deriv;
            }

            Eigen::MatrixXd derivative(const double &x, const double &y, const double &z) {
                return _derivative(x, y, z, this->param);
            }
        #endif    
 };

 #endif