#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>


template<typename G>
class YMW16 : public RegularField<G, double> {
protected:
  bool DEBUG = false;
public:
  using RegularField<G, double> :: RegularField;

  YMW16() : RegularField<G, double>() {std::cout << "YMW16 Constructor Used" << std::endl;};
  ~YMW16() {std::cout << "YMW16 Destructor Used" << std::endl;};

  SumRegularField<G, double> operator+(const RegularField<G, double>& f) {
       SumRegularField<G, double> sum(*this, f);
       return sum;
     }


// Thick disc
 double t1_n1 = 0.01132;
 double t1_h1 = 1673.;



std::vector<double> at_position(const double &x, const double &y, const double &z) const  {
  // YMW16 using a different Cartesian frame from our default one
  std::vector<double> gc_pos{y, -x, z};
  // cylindrical r
  double r_cyl{std::sqrt(gc_pos[0] * gc_pos[0] + gc_pos[1] * gc_pos[1])};
  // warp
  if (r_cyl >= r_warp) {
    double theta_warp{std::atan2(gc_pos[1], gc_pos[0])};
    gc_pos[2] -= t0_gamma_w *
                 (r_cyl - r_warp) * std::cos(theta_warp);
  }
  if (vec_length(gc_pos) > 25 ) {
    return 0.;
  } else {
    double ne{0.};
    double ne_comp[8]{0.};
    double weight_localbubble{0.};
    double weight_gum{0.};
    double weight_loop{0.};
    // longitude, in deg
    const double ec_l{
        std::atan2(gc_pos[0], r0 - gc_pos[1]) / cgs::rad};
    // call structure functions
    // since in YMW16, Fermi Bubble is not actually contributing, we ignore FB
    // for thick disk
    ne_comp[1] = thick(gc_pos[2], r_cyl);
    ne_comp[2] = thin(gc_pos[2], r_cyl);
    ne_comp[3] = spiral(gc_pos[0], gc_pos[1], gc_pos[2], r_cyl);
    ne_comp[4] = galcen(gc_pos[0], gc_pos[1], gc_pos[2]);
    ne_comp[5] = gum(gc_pos[0], gc_pos[1], gc_pos[2]);
    // localbubble boundary
    const double localbubble_boundary{110. * 0.001};
    ne_comp[6] = localbubble(gc_pos[0], gc_pos[1], gc_pos[2], ec_l,
                             localbubble_boundary);
    ne_comp[7] = nps(gc_pos[0], gc_pos[1], gc_pos[2]);
    // adding up rules
    ne_comp[0] = ne_comp[1] + std::max(ne_comp[2], ne_comp[3]);
    // distance to local bubble
    const double rlb{std::sqrt(
        std::pow(((gc_pos[1] - 8.34 ) * 0.94 - 0.34 * gc_pos[2]), 2) +
        gc_pos[0] * gc_pos[0])};
    if (rlb < localbubble_boundary) { // inside local bubble
      ne_comp[0] = rlb * ne_comp[1] +
                   std::max(ne_comp[2], ne_comp[3]);
      if (ne_comp[6] > ne_comp[0]) {
        weight_localbubble = 1;
      }
    } else { // outside local bubble
      if (ne_comp[6] > ne_comp[0] and ne_comp[6] > ne_comp[5]) {
        weight_localbubble = 1;
      }
    }
    if (ne_comp[7] > ne_comp[0]) {
      weight_loop = 1;
    }
    if (ne_comp[5] > ne_comp[0]) {
      weight_gum = 1;
    }
    // final density
    ne =
        (1 - weight_localbubble) *
            ((1 - weight_gum) * ((1 - weight_loop) * (ne_comp[0] + ne_comp[4]) +
                                 weight_loop * ne_comp[7]) +
             weight_gum * ne_comp[5]) +
        (weight_localbubble) * (ne_comp[6]);
    assert(std::isfinite(ne));
    return ne;
  }
}

 //
 class YMW16Component : public RegularThermalElectronField {
 public:
   YMW16Component() {};
   virtual ~YMW16Component() {};

   virtual double evaluate_at_pos(const std::vector<double> &pos) const = 0;

   // warp (yet not implemented)
   double r_warp = 8.4;
   double r0 = 8.3;
   double t0_gamma_w=0.14;


   double t1_ad = 2500.;
   double t1_bd = 15000.;

   double gd(const double &rr) const;

   double dgd_dad(const double &rr) const;
   double dgd_dbd(const double &rr) const;

};



class YMW16ThinDisc : public YMW16Component {
public:
 YMW16ThinDisc () {};
 virtual ~YMW16ThinDisc () {};


 // Thin disc
  double t2_a2 = 1200.;
  double t2_b2 = 4000.;
  double t2_n2 = 0.404;
  double t2_k2 = 1.54;


 double evaluate_at_pos(const std::vector<double> &pos) const;

 double dThin_dad(const std::vector<double> &pos) const;
 double dThin_dbd(const std::vector<double> &pos) const;
 double dThin_dk2(const std::vector<double> &pos) const;
 double dThin_dn2(const std::vector<double> &pos) const;
 double dThin_da2(const std::vector<double> &pos) const;
 double dThin_db2(const std::vector<double> &pos) const;

};

// ymw16 thermal electrons
class YMW16ThermalElectronField : public RegularThermalElectronField {
public:
 YMW16ThermalElectronField() {};
 virtual ~YMW16ThermalElectronField() {};

// warp
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
 double t3_rmin[5] = {3.35, 3.707, 3.56, 3.670, 8.21};
 double t3_phimin[5] = {44.4, 120.0, 218.6, 330.3, 55.1};
 double t3_tpitch[5] = {11.43, 9.84, 10.38, 10.54, 2.77};
 double t3_cpitch[5] = {11.43, 9.84, 10.38, 10.54, 2.77};
 double t3_narm[5] = {0.135, 0.129, 0.103, 0.116, 0.0057};
 double t3_warm[5] = {300., 500., 300., 500., 300.};

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



double evaluate_at_pos(const std::vector<double> &pos) const;

double thick(const double &zz, const double &rr) const;
double thin(const double &zz, const double &rr) const;
double spiral(const double &xx, const double &yy, const double &zz,
              const double &rr) const;
double galcen(const double &xx, const double &yy, const double &zz) const;
double gum(const double &xx, const double &yy, const double &zz) const;
double localbubble(const double &xx, const double &yy,
                   const double &zz, const double &ll,
                   const double &Rlb) const;
double nps(const double &xx, const double &yy, const double &zz) const;
};
