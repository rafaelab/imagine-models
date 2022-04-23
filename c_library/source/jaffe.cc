#include <cmath>
#include <vector>
#include "../headers/hamunits.h"
#include "../headers/MagneticField.h"

std::vector<double> JaffeMagneticField::evaluate_at_pos(const std::vector<double> &pos) const {
    double inner_b{0};
    if (ring) {
      inner_b = ring_amp;}
    else if (bar) {
      inner_b = bar_amp;}

    std::vector<double> bhat;
    bhat = orientation(pos);
    std::vector<double> btot{0., 0., 0.};
    double scaling = radial_scaling(pos) *
           (disk_amp * disk_scaling(pos) +
            halo_amp * halo_scaling(pos));
    for(int i=0;i<bhat.size();++i) {
        	btot[i] = bhat[i] * scaling;}
    // compress factor for each arm or for ring/bar
    std::vector<double> arm;
    arm = arm_compress(pos);
    // only inner region
    if (arm.size() == 1) {
      for(int i=0;i<bhat.size();++i) {
        btot[i] += bhat[i] * arm[0] * inner_b;}
    }
    // spiral arm region
    else {
      std::vector<double> arm_amp{arm_amp1, arm_amp2, arm_amp3, arm_amp4};
      for (decltype(arm.size()) i = 0; i < arm.size(); ++i) {
        for(int j=0;j<bhat.size();++j) {
          btot[j] += bhat[j] * arm[i] * arm_amp[i];}
      }
      
    }
    return btot;
  }

  std::vector<double> JaffeMagneticField::orientation(const std::vector<double> &pos) const {
    const double r{
        sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // cylindrical frame
    const double r_lim = ring_r;
    const double bar_lim{bar_a + 0.5 * comp_d};
    const double cos_p = std::cos(arm_pitch);
    const double sin_p = std::sin(arm_pitch); // pitch angle
    std::vector<double> tmp{0., 0., 0.};
    double quadruple{1};
    if (r < 0.5 * cgs::kpc) // forbiden region
      return tmp;
    if (pos[2] > disk_z0)
      quadruple = (1 - 2 * quadruple);
    // molecular ring
    if (ring) {
      // inside spiral arm
      if (r > r_lim) {
        tmp[0] =
            (cos_p * (pos[1] / r) - sin_p * (pos[0] / r)) * quadruple; // sin(t-p)
        tmp[1] = (-cos_p * (pos[0] / r) - sin_p * (pos[1] / r)) *
                 quadruple; //-cos(t-p)
      }
      // inside molecular ring
      else {
        tmp[0] = (1 - 2 * bss) * pos[1] / r; // sin(phi)
        tmp[1] = (2 * bss - 1) * pos[0] / r; //-cos(phi)
      }
    }
    // elliptical bar (replace molecular ring)
    else if (bar) {
      const double cos_phi = std::cos(bar_phi0);
      const double sin_phi = std::sin(bar_phi0);
      const double x = cos_phi * pos[0] - sin_phi * pos[1];
      const double y = sin_phi * pos[0] + cos_phi * pos[1];
      // inside spiral arm
      if (r > bar_lim) {
        tmp[0] =
            (cos_p * (pos[1] / r) - sin_p * (pos[0] / r)) * quadruple; // sin(t-p)
        tmp[1] = (-cos_p * (pos[0] / r) - sin_p * (pos[1] / r)) *
                 quadruple; //-cos(t-p)
      }
      // inside elliptical bar
      else {
        if (y != 0) {
          const double new_x = copysign(1, y);
          const double new_y = copysign(1, y) * (x / y) *
                                bar_b * bar_b /
                                (bar_a * bar_a);
          tmp[0] =
              (cos_phi * new_x + sin_phi * new_y) * (1 - 2 * bss);
          tmp[1] = (-sin_phi * new_x + cos_phi * new_y) *
                   (1 - 2 * bss);
          // versor
          double tmp_length = std::sqrt(tmp[0]* tmp[0] + tmp[1]* tmp[1] + tmp[2]* tmp[2]);
          if (tmp_length != 0.) {
            for(int i=0;i<tmp.size();++i) {
              tmp[i] = tmp[i]/tmp_length;}
            }
        } else {
          tmp[0] = (2 * bss - 1) * copysign(1, x) * sin_phi;
          tmp[1] = (2 * bss - 1) * copysign(1, x) * cos_phi;
        }
      }
    }
    return tmp;
  }

  double JaffeMagneticField::radial_scaling(const std::vector<double> &pos) const {
    const double r2 = pos[0] * pos[0] + pos[1] * pos[1];
    // separate into 3 parts for better view
    const double s1{
        1. - std::exp(-r2 / (r_inner * r_inner))};
    const double s2{
        std::exp(-r2 / (r_scale * r_scale))};
    const double s3{
        std::exp(-r2 * r2 /
                 (r_peak * r_peak *
                  r_peak * r_peak))};
    return s1 * (s2 + s3);
  }

  std::vector<double> JaffeMagneticField::arm_compress(const std::vector<double> &pos) const {
    const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1]) /
                      comp_r};
    const double c0{1. / comp_c - 1.};
    std::vector<double> a0 = dist2arm(pos);
    const double r_scaling{radial_scaling(pos)};
    const double z_scaling{arm_scaling(pos)};
    // for saving computing time
    const double d0_inv{(r_scaling * z_scaling) / comp_d};
    double factor{c0 * r_scaling * z_scaling};
    if (r > 1) {
      double cdrop{std::pow(r, -comp_p)};
      for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
        a0[i] = factor * cdrop *
                std::exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
      }
    } else {
      for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
        a0[i] = factor * std::exp(-a0[i] * a0[i] * d0_inv * d0_inv);
      }
    }
    return a0;
  }

  std::vector<double> JaffeMagneticField::arm_compress_dust(const std::vector<double> &pos) const {
    const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1]) /
                      comp_r};
    const double c0{1. / comp_c - 1.};
    std::vector<double> a0 = dist2arm(pos);
    const double r_scaling{radial_scaling(pos)};
    const double z_scaling{arm_scaling(pos)};
    // only difference from normal arm_compress
    const double d0_inv{(r_scaling) / comp_d};
    double factor{c0 * r_scaling * z_scaling};
    if (r > 1) {
      double cdrop{std::pow(r, -comp_p)};
      for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
        a0[i] = factor * cdrop *
                std::exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
      }
    } else {
      for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
        a0[i] = factor * std::exp(-a0[i] * a0[i] * d0_inv * d0_inv);
      }
    }
    return a0;
  }

  std::vector<double> JaffeMagneticField::dist2arm(const std::vector<double> &pos) const {
    std::vector<double> d;
    const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
    const double r_lim{ring_r};
    const double bar_lim{bar_a + 0.5 * comp_d};
    const double cos_p{std::cos(arm_pitch)};
    const double sin_p{std::sin(arm_pitch)}; // pitch angle
    const double beta_inv{-sin_p / cos_p};
    double theta{atan2(pos[1], pos[0])};
    if (theta < 0)
      theta += 2 * cgs::pi;
    // if molecular ring
    if (ring) {
      // in molecular ring, return single element vector
      if (r < r_lim) {
        d.push_back(std::fabs(ring_r - r));
      }
      // in spiral arm, return vector with arm_num elements
      else {
        // loop through arms
        std::vector<double> arm_phi{arm_phi1, arm_phi2, arm_phi3, arm_phi4};
        for (int i = 0; i < arm_num; ++i) {
          double d_ang{arm_phi[i] - theta};
          double d_rad{
              std::fabs(arm_r0 * std::exp(d_ang * beta_inv) - r)};
          double d_rad_p{
              std::fabs(arm_r0 *
                            std::exp((d_ang + 2 * cgs::pi) * beta_inv) -
                        r)};
          double d_rad_m{
              std::fabs(arm_r0 *
                            std::exp((d_ang - 2 * cgs::pi) * beta_inv) -
                        r)};
          d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
        }
      }
    }
    // if elliptical bar
    else if (bar) {
      const double cos_tmp{std::cos(bar_phi0) * pos[0] / r -
                              std::sin(bar_phi0) * pos[1] /
                                  r}; // cos(phi)cos(phi0) - sin(phi)sin(phi0)
      const double sin_tmp{std::cos(bar_phi0) * pos[1] / r +
                              std::sin(bar_phi0) * pos[0] /
                                  r}; // sin(phi)cos(phi0) + cos(phi)sin(phi0)
      // in bar, return single element vector
      if (r < bar_lim) {
        d.push_back(
            std::fabs(bar_a * bar_b /
                          sqrt(bar_a * bar_a *
                                   sin_tmp * sin_tmp +
                               bar_b * bar_b *
                                   cos_tmp * cos_tmp) -
                      r));
      }
      // in spiral arm, return vector with arm_num elements
      else {
        // loop through arms
        std::vector<double> arm_phi{arm_phi1, arm_phi2, arm_phi3, arm_phi4};
        for (int i = 0; i < arm_num; ++i) {
          double d_ang{arm_phi[i] - theta};
          double d_rad{
              std::fabs(arm_r0 * std::exp(d_ang * beta_inv) - r)};
          double d_rad_p{
              std::fabs(arm_r0 *
                            std::exp((d_ang + 2 * cgs::pi) * beta_inv) -
                        r)};
          double d_rad_m{
              std::fabs(arm_r0 *
                            std::exp((d_ang - 2 * cgs::pi) * beta_inv) -
                        r)};
          d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
        }
      }
    }
    return d;
  }

  double JaffeMagneticField::arm_scaling(const std::vector<double> &pos) const {
    return 1. / (std::cosh(pos[2] / arm_z0) *
                 std::cosh(pos[2] / arm_z0));
  }

  double JaffeMagneticField::disk_scaling(const std::vector<double> &pos) const {
    return 1. / (std::cosh(pos[2] / disk_z0) *
                 std::cosh(pos[2] / disk_z0));
  }

  double JaffeMagneticField::halo_scaling(const std::vector<double> &pos) const {
    return 1. / (std::cosh(pos[2] / halo_z0) *
                 std::cosh(pos[2] / halo_z0));
  }
