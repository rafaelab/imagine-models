#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>


template<typename G>
class JF12RandomField : public RandomField<G, std::vector<double>> {
  protected:
    bool DEBUG = false;
  public:
    using RandomField<G, std::vector<double>> :: RandomField;

    JF12RandomField() : RandomField<G, std::vector<double>>() {};
    ~JF12RandomField() {};

    const double b0_1 = b0_1;
    const double b0_2 = b0_2;
    const double b0_3 = b0_3;
    const double b0_4 = b0_4;
    const double b0_5 = b0_5;
    const double b0_6 = b0_6;
    const double b0_7 = b0_7;
    const double b0_8 = b0_8;

    const double b0_int = b0_int;
    const double z0_spiral = z0_spiral;
    const double b0_halo = b0_halo;
    const double r0_halo = r0_halo;
    const double z0_halo = z0_halo;

    const double Rmax = 20.;
    const double rho_GC = 1.;

    double spatial_profile(const double &x, const double &y, const double &z) const {

      const double r{sqrt(x * x + y * y)};
      const double rho{
          sqrt(x*x + y*y + z*z)};
      const double phi{atan2(y, z)};

      const double rc_B[8] = {
          5.1, 6.3,  7.1,  8.3,
          9.8, 11.4, 12.7, 15.5}; // neg x crossings of spiral arms
      const double inc = 11.5; // inclination, in degrees

      const double b_arms[8] = {b0_1, b0_2, b0_3, b0_4, b0_5, b0_6, b0_7, b0_8};

      double scaling_disk = 0.0;
      double scaling_halo = 0.0;

      // boundaries outside which B is zero, not sure if this works?
      if (r > Rmax || rho < rho_GC) {
        return 0.0;
      }
      if (r < 5.) {
        scaling_disk = b0_int;
      } else {
        double r_negx =
            r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi - M_PI));
        if (r_negx > rc_B[7]) {
          r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + M_PI));
        }
        if (r_negx > rc_B[7]) {
          r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + 3 * M_PI));
        }
        for (int i = 7; i >= 0; i--) {
          if (r_negx < rc_B[i] ) {
            scaling_disk = b_arms[i] * (5.) / r;
          }
        } // "region 8,7,6,..,2"
      }

      scaling_disk = scaling_disk * exp(-0.5 * z * z / (z0_spiral * z0_spiral));
      scaling_halo =
          b0_halo * exp(-std::fabs(r / r0_halo)) * exp(-std::fabs(z / z0_halo));

      return (scaling_disk * scaling_disk + scaling_halo * scaling_halo);
  }


};
