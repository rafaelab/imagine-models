#include <vector>
#include <stdexcept>

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

  virtual std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const {std::vector<double> _b(3); return _b;}

  std::vector<std::vector<std::vector<std::vector<double>>>> evaluate_grid(const std::vector<double> grid_x,
                                       const std::vector<double> grid_y,
                                       const std::vector<double> grid_z) const;
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

      std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const;
 };

