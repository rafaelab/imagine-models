#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

//Sun et al. A&A V.477 2008 ASS+RING model magnetic field

struct SunParams : Params {
    //parameter_names = {"ampx", "ampy"};
    number b_B0 = 2.;
    number b_Rsun = 8.5; 
    number b_R0 = 10.;
    number b_z0 = 1.;
    number b_Rc = 5.;
    number b_Bc = 2.;
    number b_pitch_deg = -12.;
    
    number bH_B0 = 2.;
    number bH_z0 = 1.5;
    number bH_z1a = 0.2;
    number bH_z1b = 0.4;
    number bH_R0 = 4.;
    };

class SunMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
        vector  _at_position (const double &x, const double &y, const double &z, const SunParams &p) const;
    public:
        using RegularVectorField :: RegularVectorField;

        SunParams param;

        vector at_position(const double &x, const double &y, const double &z) const {
            return _at_position(x, y, z, this->param);
        }

    #if autodiff_FOUND
    Eigen::MatrixXd _derivative(const double &x, const double &y, const double &z,  SunParams &p) {
        vector out;
        Eigen::MatrixXd deriv = ad::jacobian([&](auto _x, auto _y, auto _z, auto _p) {return this->_at_position(_x, _y, _z, _p);}, ad::wrt(p.b_B0, p.b_Bc, p.b_R0, p.b_Rc, p.b_z0, p.b_Rsun, p.bH_B0, p.bH_R0, p.bH_z0, p.bH_z1a, p.bH_z1b), ad::at(x, y, z, p), out);  
        return deriv;
    }

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z) {
        return _derivative(x, y, z, this->param);
    }
    #endif

 };
