#ifndef YMW16_H
#define YMW16_H

#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>
#include <cassert> 

#include "Field.h"
#include "RegularField.h"


struct YMWParams : Params {
  double r_warp = 8.4;
  number r0 = 8.3;
  double t0_gamma_w = 0.14;

  // Thick disc
  number t1_ad = 2500.;
  number t1_bd = 15000.;
  number t1_n1 = 0.01132;
  number t1_h1 = 1673.;

  // Thin disc
  number t2_a2 = 1200.;
  number t2_b2 = 4000.;
  number t2_n2 = 0.404;
  number t2_k2 = 1.54;

  // spiralarms
  number t3_b2s = 4000.;
  number t3_ka = 5.015;
  number t3_aa = 11680.;

  number t3_ncn = 2.4;
  number t3_wcn = 8.2;
  number t3_thetacn = 109.;
  number t3_nsg = 0.626;
  number t3_wsg = 20;
  number t3_thetasg = 78.8;
  std::array<double, 5> t3_rmin{3.35, 3.707, 3.56, 3.670, 8.21};
  std::array<double, 5> t3_phimin{44.4, 120.0, 218.6, 330.3, 55.1};
  std::array<double, 5> t3_tpitch{11.43, 9.84, 10.38, 10.54, 2.77};
  std::array<double, 5> t3_cpitch{11.43, 9.84, 10.38, 10.54, 2.77};
  std::array<double, 5> t3_narm{0.135, 0.129, 0.103, 0.116, 0.0057};
  std::array<double, 5> t3_warm{300., 500., 300., 500., 300.};

  // Galactic Center
  number t4_ngc = 6.2;
  number t4_agc = 160.;
  number t4_hgc = 35.;

  // gum
  number t5_kgn = 1.4;
  number t5_ngn = 1.84;
  number t5_wgn = 15.1;
  number t5_agn = 125.8;

  // local bubble
  number t6_j_lb = 0.480;
  number t6_nlb1 = 1.094;
  number t6_detlb1 = 28.4;
  number t6_wlb1 = 14.2;
  number t6_hlb1 = 112.9;
  number t6_thetalb1 = 195.4;
  number t6_nlb2 = 2.33;
  number t6_detlb2 = 14.7;
  number t6_wlb2 = 15.6;
  number t6_hlb2 = 43.6;
  number t6_thetalb2 = 278.2;

  // loop
  number t7_nli = 1.907;
  number t7_rli = 80.;
  number t7_wli = 15.;
  number t7_detthetali = 30.0;
  number t7_thetali = 40.0;

  };



class YMW16 : public RegularScalarField {
protected:
  bool DEBUG = false;
  number _at_position(const double &x, const double &y, const double &z, const YMWParams &p) const;
public:
  using RegularScalarField :: RegularScalarField;

  number at_position(const double &x, const double &y, const double &z) const {
        return _at_position(x, y, z, this->param);
  }

  YMWParams param;

  number thick(const double &zz, const double &rr, const YMWParams &p) const;
  number thin(const double &zz, const double &rr, const YMWParams &p) const;
  number spiral(const double &xx, const double &yy, const double &zz,
                const double &rr, const YMWParams &p) const;
  number galcen(const double &xx, const double &yy, const double &zz, const YMWParams &p) const;
  number gum(const double &xx, const double &yy, const double &zz, const YMWParams &p) const;
  number localbubble(const double &xx, const double &yy,
                    const double &zz, const double &ll,
                    const double &Rlb, const YMWParams &p) const;
  number nps(const double &xx, const double &yy, const double &zz, const YMWParams &p) const;
    
  #if autodiff_FOUND
  Eigen::VectorXd _derivative(const double &x, const double &y, const double &z,  YMWParams &p) {
      number out;
      Eigen::VectorXd deriv = ad::gradient([&](auto _x, auto _y, auto _z, auto _p) {return this->_at_position(_x, _y, _z, _p);}, ad::wrt(p.t1_h1, p.t1_n1), ad::at(x, y, z, p), out);  
      return deriv;
  }

  Eigen::VectorXd derivative(const double &x, const double &y, const double &z) {
      return _derivative(x, y, z, this->param);
  }
  #endif

};

#endif
