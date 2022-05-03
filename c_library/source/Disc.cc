
#include <cmath>


#include <hamunits.h>
#include <ThermalElectronField.h>



double YMW16ThickDisc::evaluate_at_pos(const std::vector<double> &pos) const {
    double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
    if (pos[2] > 10. * t1_h1)
      return 0.; // timesaving
    double gd = YMW16ThickDisc::gd(rr);
    return t1_n1 * gd *
           std::pow(1. / std::cosh(pos[2] / t1_h1), 2);
}



double YMW16ThickDisc::dThick_dn1(const std::vector<double> &pos) const {
    double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
    if (pos[2] > 10. * t1_h1)
      return 0.; // timesaving
    double gd = gd(rr);
    return gd * std::pow(1. / std::cosh(pos[2] / t1_h1), 2);
}

double YMW16ThickDisc::dThick_dh1(const std::vector<double> &pos) const {
    double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
    if (pos[2] > 10. * t1_h1)
      return 0.; // timesaving
    double gd = gd(rr);
    return - t1_n1 * gd *
           std::pow(1. / std::cosh(pos[2] / t1_h1), 2) *
           std::pow(std::tanh(pos[2] / t1_h1), 2) * pos[2] / t1_h1 / t1_h1;
    return evaluate;
}

double YMW16ThickDisc::dThick_dad(const std::vector<double> &pos) const {
    double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
    if (pos[2] > 10. * t1_h1)
      return 0.; // timesaving
    double dgd_dad = dgd_dad(rr);
    return dgd_dad * t1_n1 * std::pow(1. / std::cosh(pos[2] / t1_h1), 2);
}

double YMW16ThickDisc::dThick_dbd(const std::vector<double> &pos) const {
    double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
    if (pos[2] > 10. * t1_h1)
      return 0.; // timesaving
    double dgd_dbd = Disc::dgd_dbd(rr);
    return dgd_dbd * t1_n1 * std::pow(1. / std::cosh(pos[2] / t1_h1), 2);
}

////////////////////////////////////////////////////////////////////////////////


double YMW16ThinDisc::evaluate_at_pos(const std::vector<double> &pos) const {
    double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
    // z scaling, K_2*h0 in ref
    double h0{t2_k2 *
                 (32 * 0.001 + 1.6e-3 * rr + (4.e-7 / cgs::pc) * rr * rr)};
    if (pos[2] > 10. * h0)
      return 0.; // timesaving
    double gd = gd(rr);
    }
    return t2_n2 * gd *
           std::pow(1. / std::cosh((rr - t2_b2) /
                                   t2_a2),
                    2) *
           std::pow(1. / std::cosh(pos[2] / h0), 2);
}


double YMW16ThinDisc::dThin_dn2(const std::vector<double> &pos) const {
  double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  // z scaling, K_2*h0 in ref
  double h0{t2_k2 *
               (32 * 0.001 + 1.6e-3 * rr + (4.e-7 / cgs::pc) * rr * rr)};
  if (pos[2] > 10. * h0)
    return 0.; // timesaving
  double gd = gd(rr);
  }
  return gd *
         std::pow(1. / std::cosh((rr - t2_b2) /
                                 t2_a2),
                  2) *
         std::pow(1. / std::cosh(pos[2] / h0), 2);
}

double YMW16ThickDisc::dThin_dk2(const std::vector<double> &pos) const {
  double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  // z scaling, K_2*h0 in ref
  double k2_factor{32 * 0.001 + 1.6e-3 * rr + (4.e-7 / cgs::pc) * rr * rr};
  double h0{t2_k2 * k2_factor};
  if (pos[2] > 10. * h0)
    return 0.; // timesaving
  double gd = gd(rr);
  }
  return - t2_n2 * gd *
         std::pow(1. / std::cosh((rr - t2_b2) / t2_a2), 2) *
         std::pow(1. / std::cosh(pos[2] / h0), 2) *
         std::tanh(pos[2] / h0) * pos[2] / h0 / t2_k2;
}


double YMW16ThickDisc::dThick_dad(const std::vector<double> &pos) const {
  double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  // z scaling, K_2*h0 in ref
  double k2_factor{32 * 0.001 + 1.6e-3 * rr + (4.e-7 / cgs::pc) * rr * rr};
  double h0{t2_k2 * k2_factor};
  if (pos[2] > 10. * h0)
    return 0.; // timesaving
  double dgd_dad = dgd_dad(rr);
  }
  return - t2_n2 * dgd_dad *
         std::pow(1. / std::cosh((rr - t2_b2) / t2_a2), 2) *
         std::pow(1. / std::cosh(pos[2] / h0), 2) *
         std::tanh(pos[2] / h0) * pos[2] / h0 / t2_k2;
}

double YMW16ThickDisc::dThick_dbd(const std::vector<double> &pos) const {
  double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  // z scaling, K_2*h0 in ref
  double k2_factor{32 * 0.001 + 1.6e-3 * rr + (4.e-7 / cgs::pc) * rr * rr};
  double h0{t2_k2 * k2_factor};
  if (pos[2] > 10. * h0)
    return 0.; // timesaving
  double dgd_dbd = dgd_dbd(rr);
  }
  return - t2_n2 * dgd_dbd *
         std::pow(1. / std::cosh((rr - t2_b2) / t2_a2), 2) *
         std::pow(1. / std::cosh(pos[2] / h0), 2) *
         std::tanh(pos[2] / h0) * pos[2] / h0 / t2_k2;
}

double YMW16ThickDisc::dThick_db2(const std::vector<double> &pos) const {
  double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  // z scaling, K_2*h0 in ref
  double h0{t2_k2 *
               (32 * 0.001 + 1.6e-3 * rr + (4.e-7 / cgs::pc) * rr * rr)};
  if (pos[2] > 10. * h0)
    return 0.; // timesaving
  double gd = gd(rr);
  }
  return t2_n2 * gd *
         std::pow(1. / std::cosh((rr - t2_b2) /
                                 t2_a2),
                  2) *
         2 * std::tanh((rr - t2_b2) / t2_a2)) / t2_a2 *
         std::pow(1. / std::cosh(pos[2] / h0), 2);
}

double YMW16ThickDisc::dThick_da2(const std::vector<double> &pos) const {
  double rr{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  // z scaling, K_2*h0 in ref
  double h0{t2_k2 *
               (32 * 0.001 + 1.6e-3 * rr + (4.e-7 / cgs::pc) * rr * rr)};
  if (pos[2] > 10. * h0)
    return 0.; // timesaving
  double gd = gd(rr);
  }
  return t2_n2 * gd *
         std::pow(1. / std::cosh((rr - t2_b2) /
                                 t2_a2),
                  2) *
         2 * std::tanh((rr - t2_b2) / t2_a2)) / t2_a2 / t2_a2 * (rr - t2_b2)  *
         std::pow(1. / std::cosh(pos[2] / h0), 2);
}

////////////////////////////////////////////////////////////////////////////////
