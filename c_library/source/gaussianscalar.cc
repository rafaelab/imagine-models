#include <cmath>
#include <cassert>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/GaussianScalar.h"

void GaussianScalarField::_on_grid(double* val, const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) {

      fftw_complex* val_comp = construct_plans(val, shp);;
        
      draw_random_numbers_complex(val_comp, shp, inc, seed);
      
      fftw_execute(c2r);

      // normalize and add mean
      std::array<int, 3> padded_shape = {shp[0],  shp[1],  2*(shp[2]/2 + 1)}; 
      int gs_padded = grid_size(padded_shape);
      int gs = grid_size(shp);
      for (int s = 0; s < gs_padded; ++s) {
        val[s] = val[s]/std::sqrt(gs) + mean;  
      }
};


double GaussianScalarField::calculate_fourier_sigma(const double &abs_k) const {
  double sigma = simple_spectrum(abs_k, spectral_amplitude, spectral_offset, spectral_slope);
  return sigma;
}