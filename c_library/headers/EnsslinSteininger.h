#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "Field.h"
#include "RandomField.h"

class ESField : public RandomVectorField {
  protected:
    bool DEBUG = false;
  public:
    using RandomVectorField :: RandomVectorField;

    double r0 = 8.5;
    double z0 = 1.5;
    std::array<double, 3> observer{8.5, 0, 0};


    void _on_grid(std::array<double*, 3> freal,  std::array<fftw_complex*, 3> fcomp,  std::array<fftw_plan, 3> real_to_comp, std::array<fftw_plan, 3> comp_to_real, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {

      draw_random_numbers(fcomp, grid_shape, grid_increment, seed);

      for (int i =0; i<3; ++i) {
        fftw_execute(comp_to_real[i]);
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
        fftw_execute(real_to_comp[i]);
      }

      
      divergence_cleaner(fcomp[0], fcomp[1], fcomp[2], grid_shape, grid_increment);

      for (int i =0; i<3; ++i) {
        fftw_execute(comp_to_real[i]);
      }
    }


    double spatial_profile(const double &x, const double &y, const double &z) const {
      const double r_cyl{std::sqrt(x * x + y * y) -
                        std::fabs(observer[0])};
      const double zz{std::fabs(z) - std::fabs(observer[2])};
  return std::exp(-r_cyl / r0) * std::exp(-zz / z0);
}


};
