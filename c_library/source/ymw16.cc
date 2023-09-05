#include <cassert>

#include "hamunits.h"
#include "YMW.h"

double YMW16::_at_position(const double &x, const double &y, const double &z, const YMWParams &p) const  { 
  // YMW16 using a different Cartesian frame from our default one
  std::vector<double> gc_pos{y, -x, z};
  // cylindrical r
  double r_cyl{sqrt(gc_pos[0] * gc_pos[0] + gc_pos[1] * gc_pos[1])};
  // warp
  if (r_cyl >= p.r_warp) {
    double theta_warp{atan2(gc_pos[1], gc_pos[0])};
    gc_pos[2] -= p.t0_gamma_w * (r_cyl - p.r_warp) * cos(theta_warp);
  }
  double vec_length = sqrt(pow(gc_pos[0], 2) + pow(gc_pos[1], 2)+ pow(gc_pos[2], 2)); 
  if (vec_length > 25 ) {
    return 0.;
  } else {
    double ne{0.};
    double ne_comp[8]{0.};
    double weight_localbubble{0.};
    double weight_gum{0.};
    double weight_loop{0.};
    // longitude, in deg
    const double ec_l{atan2(gc_pos[0], p.r0 - gc_pos[1]) / cgs::rad};
    // call structure functions
    // since in YMW16, Fermi Bubble is not actually contributing, we ignore FB
    // for thick disk
    ne_comp[1] = thick(gc_pos[2], r_cyl, p);
    ne_comp[2] = thin(gc_pos[2], r_cyl, p);
    ne_comp[3] = spiral(gc_pos[0], gc_pos[1], gc_pos[2], r_cyl, p);
    ne_comp[4] = galcen(gc_pos[0], gc_pos[1], gc_pos[2], p);
    ne_comp[5] = gum(gc_pos[0], gc_pos[1], gc_pos[2], p);
    // localbubble boundary
    const double localbubble_boundary{110. * 0.001};
    ne_comp[6] = localbubble(gc_pos[0], gc_pos[1], gc_pos[2], ec_l,
                             localbubble_boundary, p);
    ne_comp[7] = nps(gc_pos[0], gc_pos[1], gc_pos[2], p);
    // adding up rules
    ne_comp[0] = ne_comp[1] + std::max(ne_comp[2], ne_comp[3]);
    // distance to local bubble
    const double rlb{sqrt(pow(((gc_pos[1] - 8.34 ) * 0.94 - 0.34 * gc_pos[2]), 2) + gc_pos[0] * gc_pos[0])};
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

// thick disk
double YMW16::thick(const double &zz, const double &rr, const YMWParams &p) const {
  if (zz > 10. * p.t1_h1)
    return 0.; // timesaving
  double gd{1.};
  if (rr > p.t1_bd) {
    gd = pow(1. / cosh((rr - p.t1_bd) / p.t1_ad), 2);
  }
  return p.t1_n1 * gd * pow(1. / cosh(zz / p.t1_h1), 2);
}

// thin disk
double YMW16::thin(const double &zz, const double &rr, const YMWParams &p) const {
  // z scaling, K_2*h0 in ref
  double h0{p.t2_k2 * (32 * 0.001 + 1.6e-3 * rr + (4.e-7 / 0.001) * rr * rr)};
  if (zz > 10. * h0)
    return 0.; // timesaving
  double gd{1.};
  if (rr > p.t1_bd) {
    gd = pow(1. / cosh((rr - p.t1_bd) / p.t1_ad), 2);
  }
  return p.t2_n2 * gd * pow(1. / cosh((rr - p.t2_b2) / p.t2_a2), 2) *
                        pow(1. / cosh(zz / h0), 2);
}

// spiral arms
double YMW16::spiral(const double &xx, const double &yy,
                     const double &zz, const double &rr, const YMWParams &p) const {
  // structure scaling
  double scaling{1.};
  if (rr > p.t1_bd) {
    if ((rr - p.t1_bd) > 10. * p.t1_ad)
      return 0.;
    scaling = pow( 1. / cosh((rr - p.t1_bd) / p.t1_ad), 2);
  }
  // z scaling, K_a*h0 in ref
  const double h0{p.t3_ka * (32 * 0.001 + 1.6e-3 * rr +
                             (4.e-7 / 0.001) * pow(rr, 2))};
  if (zz > 10. * h0)
    return 0.; // timesaving
  scaling *= pow(1. / cosh(zz / h0), 2);
  if ((rr - p.t3_b2s) > 10. * p.t3_aa)
    return 0.; // timesaving
  // 2nd raidus scaling
  scaling *= pow(1. / cosh((rr - p.t3_b2s) / p.t3_aa), 2);
  double smin;
  double theta{atan2(yy, xx)};
  if (theta < 0)
    theta += 2 * M_PI;
  double ne3s{0.};
  // looping through arms
  for (int i = 0; i < 4; ++i) {
    // get distance to arm center
    if (i != 4) {
      double d_phi = theta - p.t3_phimin[i];
      if (d_phi < 0) { 
        d_phi += 2. * M_PI;
      }
      double d = abs(p.t3_rmin[i] * exp(d_phi * p.t3_tpitch[i]) - rr);
      double d_p = abs(p.t3_rmin[i] * exp((d_phi + 2. * M_PI) * p.t3_tpitch[i]) - rr);
      smin = std::min(d, d_p) * p.t3_cpitch[i];
    } 
    else if (i == 4 and theta >= p.t3_phimin[i] and theta < 2) { // Local arm
      smin = abs(p.t3_rmin[i] * exp((theta + 2 * M_PI - p.t3_phimin[i]) * p.t3_tpitch[i]) - rr) * p.t3_cpitch[i];
    } else {
      continue;
    }
    if (smin > 10. * p.t3_warm[i])
      continue; // timesaving
    // accumulate density
    if (i != 2) {
      ne3s += p.t3_narm[i] * scaling * pow(1. / cosh(smin / p.t3_warm[i]), 2);
    } else if (rr > 6  and
               theta * cgs::rad > p.t3_thetacn) { // correction for Carina-Sagittarius
      const double ga =
          (1. - (p.t3_nsg) * (exp(-pow((theta * cgs::rad - p.t3_thetasg) / p.t3_wsg, 2)))) *
          (1. + p.t3_ncn) * pow(1. / cosh(smin / p.t3_warm[i]), 2);
      ne3s += p.t3_narm[i] * scaling * ga;
    } else {
      const double ga =
          (1. - (p.t3_nsg) * (exp(-pow((theta * cgs::rad - p.t3_thetasg) / p.t3_wsg, 2)))) *
          (1. + p.t3_ncn * exp(-pow((theta * cgs::rad - p.t3_thetacn) / p.t3_wcn, 2))) *
            pow(1. / cosh(smin / p.t3_warm[i]), 2);
      ne3s += p.t3_narm[i] * scaling * ga;
    }
  } // end of looping through arms
  return ne3s;
}

// galactic center
double YMW16::galcen(const double &xx, const double &yy, const double &zz, const YMWParams &p) const {
  // pos of center
  const double Xgc{50. * 0.001};
  const double Ygc{0.};
  const double Zgc{-7. * 0.001};
  const double R2gc{(xx - Xgc) * (xx - Xgc) + (yy - Ygc) * (yy - Ygc)};
  if (R2gc > 10. * p.t4_agc * p.t4_agc)
    return 0.; // timesaving
  const double Ar{exp(-R2gc / (p.t4_agc * p.t4_agc))};
  if (abs(zz - Zgc) > 10. * p.t4_hgc)
    return 0.; // timesaving
  const double Az{pow(1. / cosh((zz - Zgc) / p.t4_hgc), 2)};
  return p.t4_ngc * Ar * Az;
}

// gum nebula
double YMW16::gum(const double &xx, const double &yy, const double &zz, const YMWParams &p) const {
  if (yy < 0 or xx > 0)
    return 0.; // timesaving
  // center of Gum Nebula
  const double lc{264. * cgs::rad};
  const double bc{-4. * cgs::rad};
  const double dc{450. * 0.001};
  const double xc{dc * cos(bc) * sin(lc)};
  const double yc{p.r0 - dc * cos(bc) * cos(lc)};
  const double zc{dc * sin(bc)};
  // theta is limited in I quadrant
  const double theta{
      atan2(abs(zz - zc),
            sqrt((xx - xc) * (xx - xc) + (yy - yc) * (yy - yc)))};
  const double tantheta = tan(theta);
  // zp is positive
  double zp{(p.t5_agn * p.t5_kgn) /
            sqrt(1. + p.t5_kgn * p.t5_kgn / (tantheta * tantheta))};
  // xyp is positive
  const double xyp{zp / tantheta};
  // alpha is positive
  const double xy_dist = {
      sqrt(p.t5_agn * p.t5_agn - xyp * xyp) *
      double(p.t5_agn > xyp)};
  const double alpha{atan2(p.t5_kgn * xyp, xy_dist) +
                        theta}; // add theta, timesaving
  const double R2{(xx - xc) * (xx - xc) + (yy - yc) * (yy - yc) + (zz - zc) * (zz - zc)};
  const double r2{zp * zp + xyp * xyp};
  const double D2min{(R2 + r2 - 2. * sqrt(R2 * r2)) * sin(alpha) * sin(alpha)};
  if (D2min > 10. * p.t5_wgn * p.t5_wgn)
    return 0.;
  return p.t5_ngn * exp(-D2min / (p.t5_wgn * p.t5_wgn));
}

// local bubble
double YMW16::localbubble(const double &xx, const double &yy, const double &zz, const double &ll,
                          const double &Rlb, const YMWParams &p) const {
  if (yy < 0)
    return 0.; // timesaving
  double nel{0.};
  // r_LB in ref
  const double rLB{
      sqrt(pow(((yy - 8.34 ) * 0.94 - 0.34 * zz), 2) + pow(xx, 2))};
  // l-l_LB1 in ref
  const double dl1{
      std::min(abs(ll + 360. - p.t6_thetalb1),
               abs(p.t6_thetalb1 - (ll)))};
  if (dl1 < 10. * p.t6_detlb1 or
      (rLB - Rlb) < 10. * p.t6_wlb1 or
      zz < 10. * p.t6_hlb1) // timesaving
    nel += p.t6_nlb1 *
           pow(1. / cosh(dl1 / p.t6_detlb1), 2) *
           pow(1. / cosh((rLB - Rlb) / p.t6_wlb1), 2) *
           pow(1. / cosh(zz / p.t6_hlb1), 2);
  // l-l_LB2 in ref
  const double dl2{
      std::min(abs(ll + 360 - p.t6_thetalb2),
               abs(p.t6_thetalb2 - (ll)))};
  if (dl2 < 10. * p.t6_detlb2 or
      (rLB - Rlb) < 10. * p.t6_wlb2 or
      zz < 10. * p.t6_hlb2) // timesaving
    nel += p.t6_nlb2 *
           pow(1. / cosh(dl2 / p.t6_detlb2), 2) *
           pow(1. / cosh((rLB - Rlb) / p.t6_wlb2), 2) *
           pow(1. / cosh(zz / p.t6_hlb2), 2);
  return nel;
}

// north polar spur
double YMW16::nps(const double &xx, const double &yy, const double &zz, const YMWParams &p) const {
  if (yy < 0)
    return 0.; // timesaving
  const double theta_LI{(p.t7_thetali) * cgs::rad};
  const double x_c{-10.156 * 0.001};
  const double y_c{8106.207 * 0.001};
  const double z_c{10.467 * 0.001};
  // r_LI in ref
  const double rLI{sqrt((xx - x_c) * (xx - x_c) +
                                (yy - y_c) * (yy - y_c) +
                                (zz - z_c) * (zz - z_c))};
  const double theta{acos(((xx - x_c) * (cos(theta_LI)) +
                                   (zz - z_c) * (sin(theta_LI))) /
                                  rLI) /
                        cgs::rad};
  if (theta > 10. * p.t7_detthetali or
      (rLI - p.t7_rli) > 10. * p.t7_wli)
    return 0.; // timesaving
  return (p.t7_nli) *
         exp(-pow((rLI - p.t7_rli) / p.t7_wli, 2)) *
         exp(-pow(theta / p.t7_detthetali, 2));
}
