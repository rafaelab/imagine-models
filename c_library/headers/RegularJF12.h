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
    double b_arm_1 = 0.1;
    double b_arm_2 = 3.0;
    double b_arm_3 = -0.9;
    double b_arm_4 = -0.8;
    double b_arm_5 = -2.0;
    double b_arm_6 = -4.2;
    double b_arm_7 = 0.0;
    double b_ring = 0.1;
    double h_disk = 0.40;
    double w_disk = 0.27;
    // toroidal halo parameters
    double Bn = 1.4;
    double Bs = -1.1;
    double rn = 9.22;
    double rs = 16.7;
    double wh = 0.20;
    double z0 = 5.3;
    // X-field parameters
    double B0_X = 4.6;
    double Xtheta_const = 49;
    double rpc_X = 4.8;
    double r0_X= 2.9;
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

};

#endif
