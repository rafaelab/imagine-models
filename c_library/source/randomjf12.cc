#include <cmath>
#include <cassert>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/RandomJF12.h"


double JF12RandomField::spatial_profile(const double &x, const double &y, const double &z) const {

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

};

void JF12RandomField::_on_grid(std::array<double*, 3> grid_eval, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {
     int grid_size = grid_shape[0]*grid_shape[1]*grid_shape[2];

     std::array<fftw_complex*, 3> grid_eval_comp;
      for (int i=0; i < ndim; ++i) { 
        grid_eval_comp[i] = reinterpret_cast<fftw_complex*>(grid_eval[i]);
      }
      for (int i =0; i<3; ++i) {
        draw_random_numbers(grid_eval_comp[i], grid_shape, grid_increment, seed);
        fftw_execute(c2r[i]);
      }
      auto multiply_profile = [&](double xx, double yy, double zz) {
          int _nx = (int)((xx - grid_zeropoint[0])/grid_increment[0]); 
          int _ny = (int)((yy - grid_zeropoint[1])/grid_increment[1]);
          int _nz = (int)((zz - grid_zeropoint[2])/grid_increment[2]);
          double sp = spatial_profile(xx, yy, zz);
          int idx = _nz + grid_shape[2]*(_ny + grid_shape[1]*_nx);
          std::array<double, 3> v = {
            (grid_eval[0])[idx]*sp, 
            (grid_eval[1])[idx]*sp, 
            (grid_eval[2])[idx]*sp
            };
          return v;
        };

      evaluate_function_on_grid(grid_eval, grid_shape, grid_zeropoint, grid_increment, multiply_profile);
      for (int i =0; i<3; ++i) {
        fftw_execute(r2c[i]);
      }
      
      divergence_cleaner(grid_eval_comp[0], grid_eval_comp[1], grid_eval_comp[2], grid_shape, grid_increment);
      for (int i =0; i<3; ++i) {
        fftw_execute(c2r[i]);
        for (int s = 0; s < grid_size; ++s)
          (grid_eval[i])[s] /= grid_size;  
      }
};