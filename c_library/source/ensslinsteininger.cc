#include <cmath>
#include <cassert>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/EnsslinSteininger.h"


double ESRandomField::spatial_profile(const double &x, const double &y, const double &z) const {
      const double r_cyl{std::sqrt(x * x + y * y) -
                        std::fabs(observer[0])};
      const double zz{std::fabs(z) - std::fabs(observer[2])};
  return std::exp(-r_cyl / r0) * std::exp(-zz / z0);
};


void ESRandomField::_on_grid(std::array<double*, 3> grid_eval, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {

      int grid_size = grid_shape[0]*grid_shape[1]*grid_shape[2];

      std::array<fftw_complex*, 3> grid_eval_comp;
      for (int i=0; i < ndim; ++i) { 
        grid_eval_comp[i] = reinterpret_cast<fftw_complex*>(grid_eval[i]);
      }
        
      for (int i =0; i<3; ++i) {
        draw_random_numbers_complex(grid_eval_comp[i], grid_shape, grid_increment, seed);
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