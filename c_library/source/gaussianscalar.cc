#include <cmath>
#include <cassert>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/GaussianScalar.h"

void GaussianScalarField::_on_grid(double* grid_eval, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {

      fftw_complex* grid_eval_comp;
      grid_eval_comp = reinterpret_cast<fftw_complex*>(grid_eval);
        
      draw_random_numbers_complex(grid_eval_comp, grid_shape, grid_increment, seed);
      fftw_execute(c2r);

};