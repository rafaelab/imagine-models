#include <cmath>
#include <vector>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/Jaffe.h"

std::array<double, 3> JaffeMagneticField::_at_position(const double &x, const double &y, const double &z, const JaffeParams &p) const {
    if (x == 0. && y == 0. && z == 0.) {
        return std::array<double, 3>{0., 0., 0.};
      }
    double inner_b{0};
    if (p.ring) {
      inner_b = p.ring_amp;}
    else if (p.bar) {
      inner_b = p.bar_amp;}

    std::array<double, 3> bhat = orientation(x, y, z, p);
    std::array<double, 3> btot{0., 0., 0.};

    double scaling = radial_scaling(x, y, p) *
           (p.disk_amp * disk_scaling(z, p) +
            p.halo_amp * halo_scaling(z, p));
    
    for(int i = 0; i < bhat.size(); ++i) {
        	btot[i] = bhat[i] * scaling;}
    
    // compress factor for each arm or for ring/bar
    std::vector<double> arm = arm_compress(x, y, z, p);
    // only inner region
    if (arm.size() == 1) {
      for(int i = 0;i < bhat.size(); ++i) {
        btot[i] += bhat[i] * arm[0] * inner_b;}
    }
    // spiral arm region
    else {
      std::vector<double> arm_amp{p.arm_amp1, p.arm_amp2, p.arm_amp3, p.arm_amp4};
      for (decltype(arm.size()) i = 0; i < arm.size(); ++i) {
        for(int j = 0; j < bhat.size(); ++j) {
          btot[j] += bhat[j] * arm[i] * arm_amp[i];}
      }
    }
    return btot;
  }

  std::array<double, 3> JaffeMagneticField::orientation(const double &x, const double &y, const double &z, const JaffeParams &p) const {
    
    const double r{
        sqrt(x * x + y * y)}; // cylindrical frame
    const double r_lim = p.ring_r;
    const double bar_lim{p.bar_a + 0.5 * p.comp_d};
    const double cos_p = cos(p.arm_pitch);
    const double sin_p = sin(p.arm_pitch); // pitch angle
    
    std::array<double, 3> tmp{0., 0., 0.};
    double quadruple{1.};
    if (r < 0.5) // forbidden region
      return tmp;
    if (z > p.disk_z0)
      quadruple = (1 - 2 * quadruple);
    // molecular ring
    if (p.ring) {
      // inside spiral arm
      if (r > r_lim) {
        tmp[0] = (cos_p * (y / r) - sin_p * (x / r)) * quadruple; // sin(t-p)
        tmp[1] = (-cos_p * (x / r) - sin_p * (y / r)) * quadruple; //-cos(t-p)
      }
      // inside molecular ring
      else {
        tmp[0] = (1 - 2 * p.bss) * y / r; // sin(phi)
        tmp[1] = (2 * p.bss - 1) * x / r; //-cos(phi)
      }
    }
    // elliptical bar (replace molecular ring)
    else if (p.bar) {
      const double cos_phi = cos(p.bar_phi0);
      const double sin_phi = sin(p.bar_phi0);
      const double x = cos_phi * x - sin_phi * y;
      const double y = sin_phi * x + cos_phi * y;
      // inside spiral arm
      if (r > bar_lim) {
        tmp[0] =
            (cos_p * (y / r) - sin_p * (x / r)) * quadruple; // sin(t-p)
        tmp[1] = (-cos_p * (x / r) - sin_p * (y / r)) *
                 quadruple; //-cos(t-p)
      }
      // inside elliptical bar
      else {
        if (y != 0) {
          const double new_x = copysign(1, y);
          const double new_y = copysign(1, y) * (x / y) *
                                p.bar_b * p.bar_b /
                                (p.bar_a * p.bar_a);
          tmp[0] = (cos_phi * new_x + sin_phi * new_y) * (1 - 2 * p.bss);
          tmp[1] = (-sin_phi * new_x + cos_phi * new_y) * (1 - 2 * p.bss);
          // versor
          double tmp_length = sqrt(tmp[0]* tmp[0] + tmp[1]* tmp[1] + tmp[2]* tmp[2]);
          if (tmp_length != 0.) {
            for(int i = 0; i < tmp.size(); ++i) {
              tmp[i] = tmp[i]/tmp_length;}
            }
        } else {
          tmp[0] = (2 * p.bss - 1) * copysign(1, x) * sin_phi;
          tmp[1] = (2 * p.bss - 1) * copysign(1, x) * cos_phi;
        }
      }
    }
    return tmp;
  }

  double JaffeMagneticField::radial_scaling(const double &x, const double &y, const JaffeParams &p) const {
    const double r2 = x * x + y * y;
    // separate into 3 parts for better view
    const double s1{1. - exp(-r2 / (p.r_inner * p.r_inner))};
    const double s2{exp(-r2 / (p.r_scale * p.r_scale))};
    const double s3{exp(-r2 * r2 / (p.r_peak * p.r_peak * p.r_peak * p.r_peak))};
    return s1 * (s2 + s3);
  }

  std::vector<double> JaffeMagneticField::arm_compress(const double &x, const double &y,  const double &z, const JaffeParams &p) const {
    const double r{sqrt(x * x + y * y) / p.comp_r};
    const double c0{1. / p.comp_c - 1.};
    std::vector<double> a0 = dist2arm(x, y, p);
    const double r_scaling{radial_scaling(x, y, p)};
    const double z_scaling{arm_scaling(z, p)};
    // for saving computing time
    const double d0_inv{(r_scaling * z_scaling) / p.comp_d};
    double factor{c0 * r_scaling * z_scaling};
    if (r > 1) {
      double cdrop{pow(r, -p.comp_p)};
      for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
        a0[i] = factor * cdrop * exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
      }
    } else {
      for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
        a0[i] = factor * exp(-a0[i] * a0[i] * d0_inv * d0_inv);
      }
    }
    return a0;
  }

  std::vector<double> JaffeMagneticField::arm_compress_dust(const double &x, const double &y, const double &z, const JaffeParams &p) const {
    const double r{sqrt(x * x + y * y) / p.comp_r};
    const double c0{1. / p.comp_c - 1.};
    std::vector<double> a0 = dist2arm(x, y, p);
    const double r_scaling{radial_scaling(x, y, p)};
    const double z_scaling{arm_scaling(z, p)};
    // only difference from normal arm_compress
    const double d0_inv{(r_scaling) / p.comp_d};
    double factor{c0 * r_scaling * z_scaling};
    if (r > 1) {
      double cdrop{pow(r, -p.comp_p)};
      for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
        a0[i] = factor * cdrop * exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
      }
    } else {
      for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
        a0[i] = factor * exp(-a0[i] * a0[i] * d0_inv * d0_inv);
      }
    }
    return a0;
  }

  std::vector<double> JaffeMagneticField::dist2arm(const double &x, const double &y, const JaffeParams &p) const {
    std::vector<double> d;
    const double r{sqrt(x * x + y * y)};
    const double r_lim{p.ring_r};
    const double bar_lim{p.bar_a + 0.5 * p.comp_d};
    const double cos_p{cos(p.arm_pitch)};
    const double sin_p{sin(p.arm_pitch)}; // pitch angle
    const double beta_inv{-sin_p / cos_p};
    double theta{atan2(y, x)};
    if (theta < 0)
      theta += 2 * M_PI;
    // if molecular ring
    if (p.ring) {
      // in molecular ring, return single element vector
      if (r < r_lim) {
        d.push_back(abs(p.ring_r - r));
      }
      // in spiral arm, return vector with arm_num elements
      else {
        // loop through arms
        std::vector<double> arm_phi{p.arm_phi1, p.arm_phi2, p.arm_phi3, p.arm_phi4};
        for (int i = 0; i < p.arm_num; ++i) {
          double d_ang{arm_phi[i] - theta};
          double d_rad{
              abs(p.arm_r0 * exp(d_ang * beta_inv) - r)};
          double d_rad_p{
              abs(p.arm_r0 * exp((d_ang + 2 * M_PI) * beta_inv) - r)};
          double d_rad_m{
              abs(p.arm_r0 * exp((d_ang - 2 * M_PI) * beta_inv) -
                        r)};
          d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
        }
      }
    }
    // if elliptical bar
    else if (p.bar) {
      const double cos_tmp{cos(p.bar_phi0) * x / r - sin(p.bar_phi0) * y / r}; 
      // cos(phi)cos(phi0) - sin(phi)sin(phi0)
      const double sin_tmp{cos(p.bar_phi0) * y / r + sin(p.bar_phi0) * x / r}; 
      // sin(phi)cos(phi0) + cos(phi)sin(phi0)
      // in bar, return single element vector
      if (r < bar_lim) {
        d.push_back(abs(p.bar_a * p.bar_b / sqrt(p.bar_a * p.bar_a * sin_tmp * sin_tmp +
                        p.bar_b * p.bar_b * cos_tmp * cos_tmp) - r));
      }
      // in spiral arm, return vector with arm_num elements
      else {
        // loop through arms
        std::vector<double> arm_phi{p.arm_phi1, p.arm_phi2, p.arm_phi3, p.arm_phi4};
        for (int i = 0; i < p.arm_num; ++i) {
          double d_ang{arm_phi[i] - theta};
          double d_rad{abs(p.arm_r0 * exp(d_ang * beta_inv) - r)};
          double d_rad_p{abs(p.arm_r0 * exp((d_ang + M_PI) * beta_inv) - r)};
          double d_rad_m{abs(p.arm_r0 * exp((d_ang - 2 * M_PI) * beta_inv) -
                        r)};
          d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
        }
      }
    }
    return d;
  }

  double JaffeMagneticField::arm_scaling(const double &z, const JaffeParams &p) const {
    return 1. / (cosh(z / p.arm_z0) *
                 cosh(z / p.arm_z0));
  }

  double JaffeMagneticField::disk_scaling(const double &z, const JaffeParams &p) const {
    return 1. / (cosh(z / p.disk_z0) *
                 cosh(z / p.disk_z0));
  }

  double JaffeMagneticField::halo_scaling(const double &z, const JaffeParams &p) const {
    return 1. / (cosh(z / p.halo_z0) *
                 cosh(z / p.halo_z0));
  }
