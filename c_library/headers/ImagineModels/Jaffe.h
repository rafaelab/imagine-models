#include <cmath>
#include <vector>
#include <iostream>

#include "Field.h"
#include "RegularField.h"
#include "param.h"

struct JaffeParams : Params {
      bool quadruple = false; // quadruple pattern in halo
      bool bss = false; // bi-symmetric

      number disk_amp = 0.167;  // disk amplitude, microG
      number disk_z0 = 0.1; // disk height scale, kpc
      number halo_amp = 1.38; // halo amplitude, microG
      number halo_z0 = 3.0; // halo height scale, kpc
      number r_inner = 0.5; // inner R scale, kpc
      number r_scale = 20.; // R scale, kpc
      number r_peak = 0.; // R peak, kpc

      bool ring = false; // molecular ring
      bool bar = true; // elliptical bar
       // either ring or bar!
      number ring_amp = 0.023; // ring field amplitude, microG
      number ring_r = 5.0; // ring radius, kpc
      number bar_amp = 0.023; // bar field amplitude, microG
      number bar_a = 5.0; // major scale, kpc
      number bar_b = 3.0; // minor scale, kpc
      number bar_phi0 = 45.0; // bar major direction

      int arm_num = 4; // # of spiral arms
      number arm_r0 = 7.1; // arm ref radius, kpc
      number arm_z0 = 0.1; // arm heigth scale, kpc
      number arm_phi1 = 70; //arm ref angles, deg
      number arm_phi2 = 160;
      number arm_phi3 = 250;
      number arm_phi4 = 340;
      number arm_amp1 = 2; //  arm field amplitudes, microG
      number arm_amp2 = 0.133;
      number arm_amp3 = -3.78;
      number arm_amp4 = 0.32;
      number arm_pitch = 11.5; // pitch angle, deg

      number comp_c = 0.5; // compress factor
      number comp_d = 0.3; // arm cross-sec scale, kpc
      number comp_r = 12; // radial cutoff scale, kpc
      number comp_p = 3; //cutoff power
    };


class JaffeMagneticField : public RegularVectorField {
    protected:
      
    bool DEBUG = false;
    vector _at_position(const double &x, const double &y, const double &z, const JaffeParams &p) const;
    
    public:

    using RegularVectorField :: RegularVectorField;

    JaffeParams param;

    vector at_position(const double &x, const double &y, const double &z) const {
        return _at_position(x, y, z, this->param);
    }

    #if autodiff_FOUND
    Eigen::MatrixXd _derivative(const double &x, const double &y, const double &z,  JaffeParams &p) {
        vector out;
        Eigen::MatrixXd deriv = ad::jacobian([&](auto _x, auto _y, auto _z, auto _p) {return this->_at_position(_x, _y, _z, _p);}, ad::wrt(p.arm_amp1), ad::at(x, y, z, p), out);  
        return deriv;
    }

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z) {
        return _derivative(x, y, z, this->param);
    }
    #endif

    std::array<double, 3> at_position(const double &x, const double &y, const double &z, const JaffeParams &p) const;

    vector orientation(const double &x, const double &y, const double &z, const JaffeParams &p) const;

    number radial_scaling(const double &x, const double &y, const JaffeParams &p) const;

    vector arm_compress(const double &x, const double &y,  const double &z,const JaffeParams &p) const;

    vector arm_compress_dust(const double &x, const double &y, const double &z, const JaffeParams &p) const;

    vector dist2arm(const double &x, const double &y, const JaffeParams &p) const;

    number arm_scaling(const double &z, const JaffeParams &p) const;

    number disk_scaling(const double &z, const JaffeParams &p) const;

    number halo_scaling(const double &z, const JaffeParams &p) const;
};
