#include <cmath>
#include <cassert>
#include <iostream>
#include "ImagineModels/hamunits.h"
#include "ImagineModelsRandom/LogNormal.h"

void LogNormalScalarField::_on_grid(double* val, const std::array<int, 3> &shp, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {

      fftw_complex* val_comp = construct_plans(val, shp);;
        
      draw_random_numbers_complex(val_comp, shp, grid_increment, seed);
      
      fftw_execute(c2r);

      // normalize, add mean and exponentiate
      int gs = grid_size(shp);
      for (int s = 0; s < gs; ++s)
        val[s] = std::exp(val[s]/std::sqrt(gs) + log_mean);  
};


double LogNormalScalarField::calculate_fourier_sigma(const double &abs_k) const {
  double sigma = simple_spectrum(abs_k, spectral_amplitude, spectral_offset, spectral_slope);
  return sigma;
}