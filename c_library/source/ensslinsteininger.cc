#include <cmath>
#include <cassert>
#include <iostream>
#include "ImagineModels/hamunits.h"
#include "ImagineModelsRandom/EnsslinSteininger.h"


double ESRandomField::calculate_fourier_sigma(const double &abs_k) const {
  double sigma = simple_spectrum(abs_k, spectral_amplitude, spectral_offset, spectral_slope);
  return sigma;
};

double ESRandomField::spatial_profile(const double &x, const double &y, const double &z) const {
      const double r_cyl{std::sqrt(x * x + y * y)};
      const double zz{std::fabs(z)};
  return std::exp(-r_cyl / r0) * std::exp(-zz / z0);
};


void ESRandomField::_on_grid(std::array<double*, 3> val, const std::array<int, 3> &shp, const std::array<double, 3> &zpt, const std::array<double, 3> &inc, const int seed) {

      int gs = grid_size(shp);
      
      std::array<fftw_complex*, 3> val_comp = construct_plans(val, shp); 
      auto gen_int = std::mt19937(seed);
      std::uniform_int_distribution<int> uni(0, 1215752192);
        
      for (int i =0; i<3; ++i) {
        int sub_seed = uni(gen_int); 
        draw_random_numbers_complex(val_comp[i], shp, inc, sub_seed);
        fftw_execute(c2r[i]);
      }

      auto multiply_profile = [&](double xx, double yy, double zz) {
          int _nx = (int)((xx - zpt[0])/inc[0]); 
          int _ny = (int)((yy - zpt[1])/inc[1]);
          int _nz = (int)((zz - zpt[2])/inc[2]);
          double sp = spatial_profile(xx, yy, zz);
          int idx = _nz + shp[2]*(_ny + shp[1]*_nx);
          std::array<double, 3> v = {
            (val[0])[idx]*sp, 
            (val[1])[idx]*sp, 
            (val[2])[idx]*sp
            };
          return v;
        };
      std::array<int, 3> padded_shape = {shp[0],  shp[1],  2*(shp[2]/2 + 1)}; 
      int padded_size = padded_shape[0]*padded_shape[1]*padded_shape[2];
      evaluate_function_on_grid(val, padded_shape, zpt, inc, multiply_profile);

      for (int i =0; i<3; ++i) {
        fftw_execute(r2c[i]);
      }
      
      divergence_cleaner(val_comp[0], val_comp[1], val_comp[2], shp, inc);

      for (int i =0; i<3; ++i) {
        fftw_execute(c2r[i]);
        for (int s = 0; s < padded_size; ++s)
          (val[i])[s] /= gs;  
      }

};