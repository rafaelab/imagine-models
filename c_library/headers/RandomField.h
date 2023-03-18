#ifndef RANDOMFIELD_H
#define RANDOMFIELD_H

#include <vector>
#include <cassert>
#include <functional>
#include <iostream>
#include <algorithm>
#include <complex> 
#include <fftw3.h>
#include <random>
#include <memory>
#include <initializer_list>

#include "exceptions.h"
#include "Field.h"



template<typename POSTYPE, typename GRIDTYPE>
class RandomField : public Field<POSTYPE, GRIDTYPE>  {
protected:

    bool created_fftw_plans = false;

public:
  // Constructors
  using Field<POSTYPE, GRIDTYPE> :: Field;

  // Fields
  int seed = 0;

  // TODO: right now these parameters are not used (fixed spectrum in draw random numbers)
  double rms = 1;
  double k0 = 10; // 1/injection sale, 1/kpc
  double k1= 0.1; //  inverse cascading cutoff, 1/kpc
  double a0 = 1.7; // Kolmogorov
  double a1 = 0; // inverse cascade
  double rho = 0; // [0, inf[

  // methods
  POSTYPE at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
    // Here comes the interpolator
  }
  
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  void draw_random_numbers_complex(fftw_complex* vec,  const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_increment, const int seed); 

  GRIDTYPE on_grid(const std::vector<double>  &grid_x, const std::vector<double>  &grid_y, const std::vector<double>  &grid_z) {
    throw NotImplementedException();
  }

  double simple_spectrum(const double &abs_k, const double &rms, const double &k0, const double &s) const {
      const double p0 = rms*rms;
      double pi = 3.141592653589793;
      const double unit = 1. / (4 * pi * abs_k * abs_k);   // units fixing, wave vector in 1/kpc units
      const double P = 1 / std::pow((abs_k + k0), s);
      return P * p0 * unit;
      }

  // this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
  // original author: https://github.com/gioacchinowang
  double hammurabi_spectrum(const double &abs_k, const double &rms, const double &k0, const double &k1, const double &a0, const double &a1) const {
      const double p0 = rms*rms;
      double pi = 3.141592653589793;
      const double unit = 1. / (4 * pi * abs_k * abs_k);   // units fixing, wave vector in 1/kpc units
      // power laws
      const double band1 = double(abs_k < k1);
      const double band2 = double(abs_k > k1) * double(abs_k < k0);
      const double band3 = double(abs_k > k0);
      const double P = band1 * std::pow(k0 / k1, a1) * std::pow(abs_k / k1, 6.0) +
                      band2 / std::pow(abs_k / k0, a1) +
                      band3 / std::pow(abs_k / k0, a0);
      return P * p0 * unit;
      }

};


//include the random field method implementations (at this position, due to the use of templates)
#include "random.tpp"


class RandomScalarField : public RandomField<double, double*>  {
protected:
    // Fields
    fftw_plan c2r;
    fftw_plan r2c;
    // fftw_complex* grid_eval_comp = NULL;
    // constructors
    //RandomField() : Field<G, T>() {}
    // methods

  double* allocate_memory(std::array<int, 3> shp) {
    double* grid_eval = (double*) fftw_alloc_real(shp[0]*shp[1]*2*(shp[2]/2 + 1));
    return grid_eval;  
  }  

  void free_memory(double* grid_eval) {
      fftw_free(grid_eval);
  }

  void construct_plans(double* grid_eval, std::array<int, 3> shp) {
    fftw_complex* grid_eval_comp = reinterpret_cast<fftw_complex*>(grid_eval);
    r2c = fftw_plan_dft_r2c_3d(shp[0], shp[1], shp[2], grid_eval, grid_eval_comp, FFTW_ESTIMATE);
    c2r = fftw_plan_dft_c2r_3d(shp[0], shp[1], shp[2], grid_eval_comp, grid_eval,  FFTW_ESTIMATE);
    created_fftw_plans = true;
  }

  void destroy_plans() {
    if (created_fftw_plans) {
      fftw_destroy_plan(c2r);
      fftw_destroy_plan(r2c);
    }

  }

public:
  RandomScalarField() : RandomField<double, double*>() {
    };

  RandomScalarField(std::array<int, 3>  shape, std::array<double, 3>  zeropoint, std::array<double, 3>  increment) : RandomField<double, double*>(shape, zeropoint, increment) {
    //accumulate wisdom
    double* grid_eval = allocate_memory(shape);
    fftw_complex* grid_eval_comp = reinterpret_cast<fftw_complex*>(grid_eval);
    fftw_plan r2c_temp = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], grid_eval, grid_eval_comp, FFTW_ESTIMATE);
    fftw_plan c2r_temp = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], grid_eval_comp, grid_eval,  FFTW_ESTIMATE);
    fftw_destroy_plan(c2r_temp);
    fftw_destroy_plan(r2c_temp); //delete plans, but wisdom will be kept!
    fftw_free(grid_eval);
    };

  ~RandomScalarField() {
    //delete wisdom
    destroy_plans();
    #ifdef _OPENMP
          fftw_cleanup_threads();
    #else
          fftw_cleanup();
    #endif
    };

  // Fields

  const int ndim = 1;
  // methods

  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  double* on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, int seed = 0) {
    throw NotImplementedException();
    //std::vector<double> c{0};
    //return c;
  }

  // This function is the place where the global routine should be implemented, i.e. how the spatial  profile is connected to the random number, and if divergence cleaning needs to be performed. 
  virtual void _on_grid(double* grid_eval,  const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) = 0;

  double* on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {
    if (initialized_with_grid) 
      throw GridException();
    double* grid_eval = allocate_memory(grid_shape);
    construct_plans(grid_eval, grid_shape);
    _on_grid(grid_eval, grid_shape, grid_zeropoint, grid_increment, seed);
    return grid_eval;
  }

  double* on_grid(const int seed) {
    if (not initialized_with_grid) 
      throw GridException();
    double* grid_eval = allocate_memory(shape);
    construct_plans(grid_eval, shape);
    _on_grid(grid_eval,  shape, zeropoint, increment, seed);
    return grid_eval;
  }


};



class RandomVectorField : public RandomField<std::array<double, 3>, std::array<double*, 3>>  {
protected:
    // Fields

    std::array<fftw_plan, 3> r2c;
    std::array<fftw_plan, 3> c2r;

    // constructors
    //RandomField() : Field<G, T>() {}
    // methods


  std::array<double*, 3> allocate_memory(std::array<int, 3> shp) {
    std::array<double*, 3> grid_eval;
      for (int i=0; i < ndim; ++i) {
        grid_eval[i] = (double*) fftw_alloc_real(shp[0]*shp[1]*2*(shp[2]/2 + 1));
        }
    return grid_eval;  
  }  

  void free_memory(std::array<double*, 3> grid_eval) {
    for (int i=0; i < ndim; ++i) { 
      fftw_free(grid_eval[i]);
    }
  }

  void construct_plans(std::array<double*, 3> grid_eval, std::array<int, 3> shp) {
    std::array<fftw_complex*, 3> grid_eval_comp;
      for (int i=0; i < ndim; ++i) {
        grid_eval_comp[i] = reinterpret_cast<fftw_complex*>(grid_eval[i]);
        r2c[i] = fftw_plan_dft_r2c_3d(shp[0], shp[1], shp[2], grid_eval[i], grid_eval_comp[i], FFTW_ESTIMATE);
        c2r[i] = fftw_plan_dft_c2r_3d(shp[0], shp[1], shp[2], grid_eval_comp[i], grid_eval[i],  FFTW_ESTIMATE);
      }
    created_fftw_plans = true;
  }

  void destroy_plans() {
    if (created_fftw_plans) {
      for (int i=0; i < ndim; ++i) {
        fftw_destroy_plan(c2r[i]);
        fftw_destroy_plan(r2c[i]);
      }
    }
  }

public:


  RandomVectorField() : RandomField() {
    };
  RandomVectorField(std::array<int, 3>  shape, std::array<double, 3>  zeropoint, std::array<double, 3>  increment) : RandomField(shape, zeropoint, increment) {
    double* grid_eval = (double*) fftw_alloc_real(shape[0]*shape[1]*2*(shape[2]/2 + 1));
    fftw_complex* grid_eval_comp = reinterpret_cast<fftw_complex*>(grid_eval);
    fftw_plan r2c_temp = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], grid_eval, grid_eval_comp, FFTW_MEASURE);
    fftw_plan c2r_temp = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], grid_eval_comp, grid_eval,  FFTW_MEASURE);
    fftw_free(grid_eval);
    fftw_destroy_plan(c2r_temp);
    fftw_destroy_plan(r2c_temp); // plains are destroyd but wisdom is saved!
    };
  ~RandomVectorField() {
    destroy_plans();
    };

  // fields
  const int ndim = 3;
  // methods
  std::array<double, 3> at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
   // T c{0};
   // return c;
  }
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  // This function is the place where the global routine should be implemented, i.e. how the spatial  profile is connected to the random number, and if divergence cleaning needs to be performed. 
  virtual void _on_grid(std::array<double*, 3> grid_eval, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) = 0;

  std::array<double*, 3> on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, int seed = 0) {
    throw NotImplementedException();
    //std::vector<double> c{0};
    //return c;
  }

  std::array<double*, 3> on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, int nseed) {
    if (initialized_with_grid) 
      throw GridException();
    std::array<double*, 3> grid_eval = allocate_memory(grid_shape);
    construct_plans(grid_eval, grid_shape);
    _on_grid(grid_eval, grid_shape, grid_zeropoint, grid_increment, seed);
    return grid_eval;
    
  }

  std::array<double*, 3> on_grid(int nseed) {
    if (not initialized_with_grid) 
      throw GridException();
    std::array<double*, 3> grid_eval = allocate_memory(shape);
    construct_plans(grid_eval, shape);
    _on_grid(grid_eval,  shape, zeropoint, increment, nseed);
    return grid_eval;
  }

  
// this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
// original author: https://github.com/gioacchinowang
  void divergence_cleaner(fftw_complex* bx, fftw_complex* by, fftw_complex* bz,  const std::array<int, 3> &shape, const std::array<double, 3> &increment) const {
    #ifdef _OPENMP
      #pragma omp parallel for schedule(static)
    #endif
      for (int i = 0; i < shape[0]; ++i) {
        double kx = i / increment[0];
        if (i >= (shape[0] + 1) / 2)
          kx -= 1. /  increment[0];
          // it's faster to calculate indices manually
        const int idx_lv1 = i * shape[1] * shape[2];
        for (int j = 0; j < shape[1]; ++j) {
          double ky = j /  increment[1];
          if (j >= (shape[1] + 1) / 2)
            ky -= 1. /  increment[1];
          const int idx_lv2 = idx_lv1 + j * shape[2];
          for (int l = 0; l < shape[2]/2 + 1; ++l) {
            // 0th term is fixed to zero in allocation
            if (i == 0 and j == 0 and l == 0)
              continue;
            double kz = 1. /  increment[2];
            const int idx = idx_lv2 + l;
            double k_length = 0;
            double b_length = 0;
            double b_dot_k = 0;
            std::array<double, 3> k{kx, ky, kz};
            std::array<double, 3> b{(*bx)[idx], (*by)[idx], (*bz)[idx]};
            b_length = static_cast<double>(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
            k_length = static_cast<double>(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
            b_dot_k = (b[0]*k[0] + b[1]*k[1] + b[2]*k[2]);

            if (k_length == 0 or b_length == 0) {
              continue;
              }
            k_length = std::sqrt(k_length);

            const double bk_over_k = b_dot_k / k_length;
            // multiply \sqrt(3) for preserving spectral power statistically
            (*bx)[idx] = 1.73205081 * ((*bx)[idx] - k[0] * bk_over_k);
            (*by)[idx] = 1.73205081 * ((*by)[idx] - k[1] * bk_over_k);
            (*bz)[idx] = 1.73205081 * ((*bz)[idx] - k[2] * bk_over_k);
          } // l
        } // j
      } // i
    }
  };


#endif /* RANDOMFIELD_H */
