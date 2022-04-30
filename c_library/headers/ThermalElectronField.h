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


 // ymw16 thermal electrons
 class YMW16ThermalElectronField : public ThermalElectronField {
 public:
   YMW16ThermalElectronField() = default;
   virtual ~YMW16ThermalElectronField() = default;

   double r_warp, r0;
   double t0_gamma_w;
   double t1_ad, t1_bd, t1_n1, t1_h1;
   double t2_a2, t2_b2, t2_n2, t2_k2;
   double t3_b2s, t3_ka, t3_narm[5], t3_warm[5], t3_aa, t3_ncn, t3_wcn,
      t3_thetacn, t3_nsg, t3_wsg, t3_thetasg, t3_rmin[5], t3_phimin[5],
      t3_tpitch[5], t3_cpitch[5];
   double t4_ngc, t4_agc, t4_hgc;
   double t5_kgn, t5_ngn, t5_wgn, t5_agn;
   double t6_j_lb, t6_nlb1, t6_detlb1, t6_wlb1, t6_hlb1, t6_thetalb1,
      t6_nlb2, t6_detlb2, t6_wlb2, t6_hlb2, t6_thetalb2;
   double t7_nli, t7_rli, t7_wli, t7_detthetali, t7_thetali;
