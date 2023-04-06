#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>
#include <cassert> 

#include "Field.h"
#include "RegularField.h"

struct YMWParams : Params {
    //parameter_names = {"ampx", "ampy"};

  double r_warp = 8.4;
  double r0 = 8.3;
  double t0_gamma_w=0.14;

  // Thick disc
  double t1_ad = 2500.;
  double t1_bd = 15000.;
  double t1_n1 = 0.01132;
  double t1_h1 = 1673.;

  // Thin disc
  double t2_a2 = 1200.;
  double t2_b2 = 4000.;
  double t2_n2 = 0.404;
  double t2_k2 = 1.54;

  // spiralarms
  double t3_b2s = 4000.;
  double t3_ka = 5.015;
  double t3_aa = 11680.;

  double t3_ncn = 2.4;
  double t3_wcn = 8.2;
  double t3_thetacn = 109.;
  double t3_nsg = 0.626;
  double t3_wsg = 20;
  double t3_thetasg = 78.8;
  std::array<double, 5> t3_rmin{3.35, 3.707, 3.56, 3.670, 8.21};
  std::array<double, 5> t3_phimin{44.4, 120.0, 218.6, 330.3, 55.1};
  std::array<double, 5> t3_tpitch{11.43, 9.84, 10.38, 10.54, 2.77};
  std::array<double, 5> t3_cpitch{11.43, 9.84, 10.38, 10.54, 2.77};
  std::array<double, 5> t3_narm{0.135, 0.129, 0.103, 0.116, 0.0057};
  std::array<double, 5> t3_warm{300., 500., 300., 500., 300.};

  // Galactic Center
  double t4_ngc = 6.2;
  double t4_agc = 160.;
  double t4_hgc = 35.;

  // gum
  double t5_kgn = 1.4;
  double t5_ngn = 1.84;
  double t5_wgn = 15.1;
  double t5_agn = 125.8;

  // local bubble
  double t6_j_lb = 0.480;
  double t6_nlb1 = 1.094;
  double t6_detlb1 = 28.4;
  double t6_wlb1 = 14.2;
  double t6_hlb1 = 112.9;
  double t6_thetalb1 = 195.4;
  double t6_nlb2 = 2.33;
  double t6_detlb2 = 14.7;
  double t6_wlb2 = 15.6;
  double t6_hlb2 = 43.6;
  double t6_thetalb2 = 278.2;

  // loop
  double t7_nli = 1.907;
  double t7_rli = 80.;
  double t7_wli = 15.;
  double t7_detthetali = 30.0;
  double t7_thetali = 40.0;
    };

class YMW16 : public RegularScalarField {
protected:
  bool DEBUG = false;
  double _at_position(const double &x, const double &y, const double &z, const YMWParams &p) const;

public:
  using RegularScalarField :: RegularScalarField;

  YMWParams param;

  double thick(const double &zz, const double &rr, const YMWParams &p) const;
  double thin(const double &zz, const double &rr, const YMWParams &p) const;
  double spiral(const double &xx, const double &yy, const double &zz,
                const double &rr, const YMWParams &p) const;
  double galcen(const double &xx, const double &yy, const double &zz, const YMWParams &p) const;
  double gum(const double &xx, const double &yy, const double &zz, const YMWParams &p) const;
  double localbubble(const double &xx, const double &yy,
                    const double &zz, const double &ll,
                    const double &Rlb, const YMWParams &p) const;
  double nps(const double &xx, const double &yy, const double &zz, const YMWParams &p) const;

  double at_position(const double &x, const double &y, const double &z) const {
        return _at_position(x, y, z, this->param);
        }

  
};
