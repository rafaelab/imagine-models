#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "Field.h"
#include "RandomField.h"

class JF12RandomField : public RandomVectorField {
  protected:
    bool DEBUG = false;
  public:
    using RandomVectorField :: RandomVectorField;

    JF12RandomField() : RandomVectorField() {};
    ~JF12RandomField() {};

    double b0_1 = 10.81; // uG
    double b0_2 = 6.96; // uG
    double b0_3 = 9.59; // uG
    double b0_4 = 6.96; // uG
    double b0_5 = 1.96; // uG
    double b0_6 = 16.34; // uG
    double b0_7 = 37.29; // uG
    double b0_8 = 10.35; // uG
    double b0_int = 7.63; // uG
    double z0_spiral = 0.61; // kpc

    double b0_halo = 4.68; // uG
    double r0_halo = 10.97; // kpc
    double z0_halo = 2.84; // kpc

    double Rmax = 20.;
    double rho_GC = 1.;

    void _on_grid(std::array<double*, 3> freal,  std::array<fftw_complex*, 3> fcomp,  std::array<fftw_plan, 3> forward, std::array<fftw_plan, 3> backward, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {

      draw_random_numbers(fcomp, grid_shape, grid_increment, seed);

      for (int i =0; i<3; ++i) {
        fftw_execute(c2r[i]);
      }

      auto multiply_profile = [&](double xx, double yy, double zz) {
          int _nx = (int)((xx - grid_zeropoint[0])/grid_increment[0]); 
          int _ny = (int)((yy - grid_zeropoint[1])/grid_increment[1]);
          int _nz = (int)((zz - grid_zeropoint[2])/grid_increment[2]);
          double sp = spatial_profile(xx, yy, zz);
          int idx = _nz + grid_shape[2]*(_ny + grid_shape[1]*_nx);
          std::array<double, 3> v = {
            (freal[0])[idx]*sp, 
            (freal[1])[idx]*sp, 
            (freal[2])[idx]*sp
            };
          return v;
        };

      evaluate_function_on_grid(freal, grid_shape, grid_zeropoint, grid_increment, multiply_profile);

      for (int i =0; i<3; ++i) {
        fftw_execute(r2c[i]);
      }

      
      divergence_cleaner(fcomp[0], fcomp[1], fcomp[2], grid_shape, grid_increment);

      for (int i =0; i<3; ++i) {
        fftw_execute(c2r[i]);
      }
    }


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
