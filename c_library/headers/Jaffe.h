#include <cmath>
#include <vector>
#include <iostream>

#include "Field.h"
#include "RegularField.h"

class JaffeMagneticField : public RegularVectorField {
    protected:
      bool DEBUG = false;
    public:

    using RegularVectorField :: RegularVectorField;

      bool quadruple = false; // quadruple pattern in halo
      bool bss = false; // bi-symmetric

      double disk_amp = 0.167;  // disk amplitude, microG
      double disk_z0 = 0.1; // disk height scale, kpc
      double halo_amp = 1.38; // halo amplitude, microG
      double halo_z0 = 3.0; // halo height scale, kpc
      double r_inner = 0.5; // inner R scale, kpc
      double r_scale = 20.; // R scale, kpc
      double r_peak = 0.; // R peak, kpc

      bool ring = false; // molecular ring
      bool bar = true; // elliptical bar
       // either ring or bar!
      double ring_amp = 0.023; // ring field amplitude, microG
      double ring_r = 5.0; // ring radius, kpc
      double bar_amp = 0.023; // bar field amplitude, microG
      double bar_a = 5.0; // major scale, kpc
      double bar_b = 3.0; // minor scale, kpc
      double bar_phi0 = 45.0; // bar major direction

      int arm_num = 4; // # of spiral arms
      double arm_r0 = 7.1; // arm ref radius, kpc
      double arm_z0 = 0.1; // arm heigth scale, kpc
      double arm_phi1 = 70; //arm ref angles, deg
      double arm_phi2 = 160;
      double arm_phi3 = 250;
      double arm_phi4 = 340;
      double arm_amp1 = 2; //  arm field amplitudes, microG
      double arm_amp2 = 0.133;
      double arm_amp3 = -3.78;
      double arm_amp4 = 0.32;
      double arm_pitch = 11.5; // pitch angle, deg

      double comp_c = 0.5; // compress factor
      double comp_d = 0.3; // arm cross-sec scale, kpc
      double comp_r = 12; // radial cutoff scale, kpc
      double comp_p = 3; //cutoff power

    std::array<double, 3> at_position(const double &x, const double &y, const double &z) const;

    std::array<double, 3> orientation(const double &x, const double &y, const double &z) const;

    double radial_scaling(const double &x, const double &y) const;

    std::vector<double> arm_compress(const double &x, const double &y,  const double &z) const;

    std::vector<double> arm_compress_dust(const double &x, const double &y, const double &z) const;

    std::vector<double> dist2arm(const double &x, const double &y) const;

    double arm_scaling(const double &z) const;

    double disk_scaling(const double &z) const;

    double halo_scaling(const double &z) const;
/*
      std::vector<double> orientation(const double &x, const double &y,  const double &z) const;
      double radial_scaling(const double &x, const double &y) const;
      double arm_scaling(const double &z) const;
      double disk_scaling(const double &z) const;
      double halo_scaling(const double &z) const;
      std::vector<double> arm_compress(const double &x, const double &y,  const double &z) const;
      std::vector<double> arm_compress_dust(const double &x, const double &y,  const double &z) const;
      std::vector<double> dist2arm(const double &x, const double &y) const;
*/
};
