#include <vector>
#include <stdexcept>
#include <functional>

// thermal electron base class
class ThermalElectronField {
public:
  ThermalElectronField() = default;
  virtual ~ThermalElectronField() = default;
};

// regular thermal electrons
class RegularThermalElectronField : public ThermalElectronField {
public:
  RegularThermalElectronField() = default;
  virtual ~RegularThermalElectronField() = default;

  virtual double evaluate_at_pos(const std::vector<double> &pos) const = 0;

  std::vector<double> _evaluate_grid(const std::vector<double> grid_x,
                                     const std::vector<double> grid_y,
                                     const std::vector<double> grid_z,
                                     std::function<double(std::vector<double>)> ev_at_pos) const;

  std::vector<double> evaluate_grid(const std::vector<double> grid_x,
                                    const std::vector<double> grid_y,
                                    const std::vector<double> grid_z) const {
                                    return _evaluate_grid(grid_x, grid_y, grid_z, [this](std::vector<double> p) {return evaluate_at_pos(p);});
                                    };


 //
 class YMW16Component : public ThermalElectronField {
 public:
   YMW16Component() = default;
   virtual ~YMW16Component() = default;

   // warp (yet not implemented)
   double r_warp = 8.4;
   double r0 = 8.3;
   double t0_gamma_w=0.14;


   double t1_ad = 2500.;
   double t1_bd = 15000.;

   double gd(const double rr) const;

   double dgd_dad(const double rr) const;
   double dgd_dbd(const double rr) const;


class YMW16ThickDisc : public YMW16Component {
public:
 YMW16ThickDisc () = default;
 virtual ~YMW16ThickDisc () = default;


// Thick disc
  double t1_n1 = 0.01132;
  double t1_h1 = 1673.;


 double evaluate_at_pos(const std::vector<double> &pos) const;

 double dThick_dad(const std::vector<double> &pos) const;
 double dThick_dbd(const std::vector<double> &pos) const;
 double dThick_dh1(const std::vector<double> &pos) const;
 double dThick_dn1(const std::vector<double> &pos) const;

class YMW16ThinDisc : public YMW16Component {
public:
 YMW16ThinDisc () = default;
 virtual ~YMW16ThinDisc () = default;


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

// ymw16 thermal electrons
class YMW16ThermalElectronField : public ThermalElectronField {
public:
 YMW16ThermalElectronField() = default;
 virtual ~YMW16ThermalElectronField() = default;

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
 dpuble t3_thetasg = 78.8;
 double t3_rmin[5];
 double t3_phimin[5];
 double t3_tpitch[5];
 double t3_cpitch[5];
 double t3_narm[5];
 double t3_warm[5];
 
// Galactic Center
 double t4_ngc, t4_agc, t4_hgc;
 double t5_kgn, t5_ngn, t5_wgn, t5_agn;
 double t6_j_lb, t6_nlb1, t6_detlb1, t6_wlb1, t6_hlb1, t6_thetalb1,
    t6_nlb2, t6_detlb2, t6_wlb2, t6_hlb2, t6_thetalb2;
 double t7_nli, t7_rli, t7_wli, t7_detthetali, t7_thetali;
