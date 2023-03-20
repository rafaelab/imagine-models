#include <cmath>
#include <cassert>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/GaussianScalar.h"

void GaussianScalarField::_on_grid(double* val, const std::array<int, 3> &shp, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {

      fftw_complex* val_comp = construct_plans(val, shp);;
        
      draw_random_numbers_complex(val_comp, shp, grid_increment, seed);
      
      fftw_execute(c2r);

      // normalize and add mean
      int gs = grid_size(shp);
      for (int s = 0; s < gs; ++s)
        val[s] = val[s]/std::sqrt(gs) + mean;  
};


double GaussianScalarField::calculate_fourier_sigma(const double &abs_k) const {
  double sigma = simple_spectrum(abs_k, spectral_amplitude, spectral_offset, spectral_slope);
  return sigma;
}