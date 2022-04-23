#include <vector>
#include <stdexcept>
#include <functional>

// magnetic field base class
class MagneticField {
public:
  MagneticField() = default;
  virtual ~MagneticField() = default;
};

// regular magnetic field
class RegularMagneticField : public MagneticField {
public:
  RegularMagneticField() = default;
  virtual ~RegularMagneticField() = default;

  virtual std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const = 0;

  std::vector<std::vector<std::vector<std::vector<double>>>> _evaluate_grid(const std::vector<double> grid_x,
                                       const std::vector<double> grid_y,
                                       const std::vector<double> grid_z,
                                       std::function<std::vector<double>(std::vector<double> )> ev_at_pos) const;

  std::vector<std::vector<std::vector<std::vector<double>>>> evaluate_grid(const std::vector<double> grid_x,
                                    const std::vector<double> grid_y,
                                    const std::vector<double> grid_z) const {
                                    return _evaluate_grid(grid_x, grid_y, grid_z, [this](std::vector<double> p) {return evaluate_at_pos(p);});
                                    };


 };

 class JF12MagneticField : public RegularMagneticField {
    public:
      JF12MagneticField() = default;
      virtual ~JF12MagneticField() = default;

      double b_arm_1 = 0.1;
      double b_arm_2 = 3.0;
      double b_arm_3 = -0.9;
      double b_arm_4 = -0.8;
      double b_arm_5 = -2.0;
      double b_arm_6 = -4.2;
      double b_arm_7 = 0.0;
      double b_ring = 0.1;
      double h_disk = 0.40;
      double w_disk = 0.27;

                          // toroidal halo parameters
      double Bn = 1.4;
      double Bs = -1.1;
      double rn = 9.22;
      double rs = 16.7;
      double wh = 0.20;
      double z0 = 5.3;
      // X-field parameters
      double B0_X = 4.6;
      double Xtheta_const = 49;
      double rpc_X = 4.8;
      double r0_X= 2.9;

      std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override;
 };


 class HelixMagneticField : public RegularMagneticField {
    public:
      HelixMagneticField() = default;
      virtual ~HelixMagneticField() = default;

      double bx = 0.;
      double by = 0.;
      double bz = 0.;
      double rmax = 3.;
      double rmin = 0.;

      std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override;

      std::vector<double> _dbx_at_pos(const std::vector<double> &pos) const;
      std::vector<double> _dby_at_pos(const std::vector<double> &pos) const;
      std::vector<double> _dbz_at_pos(const std::vector<double> &pos) const;

      std::vector<std::vector<std::vector<std::vector<double>>>> dbx_grid(const std::vector<double> grid_x,
                                        const std::vector<double> grid_y,
                                        const std::vector<double> grid_z) const {
                                        return _evaluate_grid(grid_x, grid_y, grid_z, [this](std::vector<double> p)
                                        {return _dbx_at_pos(p);});
                                        };
      std::vector<std::vector<std::vector<std::vector<double>>>> dby_grid(const std::vector<double> grid_x,
                                        const std::vector<double> grid_y,
                                        const std::vector<double> grid_z) const {
                                        return _evaluate_grid(grid_x, grid_y, grid_z, [this](std::vector<double> p)
                                        {return _dby_at_pos(p);});
                                      };
      std::vector<std::vector<std::vector<std::vector<double>>>> dbz_grid(const std::vector<double> grid_x,
                                        const std::vector<double> grid_y,
                                        const std::vector<double> grid_z) const {
                                        return _evaluate_grid(grid_x, grid_y, grid_z, [this](std::vector<double> p)
                                        {return _dbz_at_pos(p);});
                                                                        };

 };

 class JaffeMagneticField : public RegularMagneticField {
    public:
      JaffeMagneticField() = default;
      virtual ~JaffeMagneticField() = default;

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
      bool bar = false; // elliptical bar
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

      std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override;

      std::vector<double> orientation(const std::vector<double> &pos) const;
      double radial_scaling(const std::vector<double> &pos) const;
      double arm_scaling(const std::vector<double> &pos) const;
      double disk_scaling(const std::vector<double> &pos) const;
      double halo_scaling(const std::vector<double> &pos) const;
      std::vector<double> arm_compress(const std::vector<double> &pos) const;
      std::vector<double> arm_compress_dust(const std::vector<double> &pos) const;
      std::vector<double> dist2arm(const std::vector<double> &pos) const;
};
