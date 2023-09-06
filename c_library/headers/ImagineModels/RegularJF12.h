#ifndef REGULARJF12_H
#define REGULARJF12_H

#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "Field.h"
#include "RegularField.h"
#include "param.h"


struct JF12RegularParams : Params {
    number b_arm_1 = 0.1;
    number b_arm_2 = 3.0;
    number b_arm_3 = -0.9;
    number b_arm_4 = -0.8;
    number b_arm_5 = -2.0;
    number b_arm_6 = -4.2;
    number b_arm_7 = 0.0;
    number b_ring = 0.1;
    number h_disk = 0.40;
    number w_disk = 0.27;
    // toroidal halo parameters
    number Bn = 1.4;
    number Bs = -1.1;
    number rn = 9.22;
    number rs = 16.7;
    number wh = 0.20;
    number z0 = 5.3;
    // X-field parameters
    number B0_X = 4.6;
    number Xtheta_const = 49;
    number rpc_X = 4.8;
    number r0_X= 2.9;
    };


class JF12MagneticField : public RegularVectorField {
  protected:
    bool DEBUG = false;

    vector _at_position(const double &x, const double &y, const double &z, const JF12RegularParams &p) const;
  
  public:
    using RegularVectorField :: RegularVectorField;

    JF12RegularParams param;

    vector at_position(const double &x, const double &y, const double &z) const {
        return _at_position(x, y, z, this->param);
  }

  #if autodiff_FOUND
    Eigen::MatrixXd _derivative(const double &x, const double &y, const double &z,  JF12RegularParams &p) {
        vector out;
        Eigen::MatrixXd deriv = ad::jacobian([&](auto _x, auto _y, auto _z, auto _p) {return this->_at_position(_x, _y, _z, _p);}, ad::wrt(p.b_arm_1, p.b_arm_2, p.b_arm_3, p.b_arm_4, p.b_arm_5, p.b_arm_6, p.b_arm_7,
        p.b_ring, p.h_disk, p.w_disk, 
        p.Bn, p.Bs, p.rn, p.rs, p.wh, p.z0,
        p.B0_X, p.Xtheta_const, p.rpc_X, p.r0_X), ad::at(x, y, z, p), out);  
        return deriv;
    }

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z) {
        return _derivative(x, y, z, this->param);
    }
  #endif

};

#endif
