#ifndef RANDOMFIELD_H
#define RANDOMFIELD_H

#include <vector>
#include <cassert>
#include <stdexcept>
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



template<typename T>
class RandomField : public Field<T>  {
protected:
    // Fields
    int ndim ;
    // Constructors
    using Field<T> :: Field;
    // RandomField() : Field<T>() {};
    //RandomField() : Field<G, T>() {}
    // methods

public:
  // Constructors
  RandomField(int shape_x, int shape_y, int shape_z) : Field<T>(), shape_x(shape_x), shape_y(shape_y), shape_z(shape_z) {
        }
  ~RandomField() {};

  // Fields
  int seed;
  // for destructor
  bool clean_switch = false;
  // methods
  T at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
    // Here comes the interpolator
    T c{0};
    return c;
  }
  
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  double on_grid(const std::vector<double>  &grid_x, const std::vector<double>  &grid_y, const std::vector<double>  &grid_z) {
    throw NotImplementedException();
    std::vector<double> c{0};
    return c;
  }

  virtual double on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) = 0;

  void draw_random_numbers(std::vector<fftw_complex*> c_fields, const std::vector<int> &n, const std::vector<double>  &increment, const int seed);

  double spectrum(const double &abs_k, const double &rms, const double &k0, const double &k1, const double &a0, const double &a1) const;

};

//include the random field method implementations (at this position, due to the use of templates)
#include "RandomField.tpp"

class RandomScalarField : public RandomField<double>  {
protected:
    // Fields
    int ndim = 1;
    int shape_x = 0;
    int shape_y = 0;
    int shape_z = 0;
    bool recreate_fftw_plans;
    fftw_plan c2r;
    fftw_plan r2c;
    double* field_real = 0;
    fftw_complex* field_comp = 0;
    // constructors
    //RandomField() : Field<G, T>() {}
    // methods

public:
  RandomScalarField() : RandomField<double>() {
    recreate_fftw_plans = true;
    };
  RandomScalarField(int shape_x, int shape_y, int shape_z) : RandomField<double>(), shape_x(shape_x), shape_y(shape_y), shape_z(shape_z) {
    recreate_fftw_plans = false;
    field_real = (double*) fftw_malloc(sizeof(fftw_complex) * shape_x*shape_y*(shape_z+2));
    field_comp = reinterpret_cast<fftw_complex*>(field_real);
    r2c = fftw_plan_dft_r2c_3d(shape_x, shape_y, shape_z, field_real, field_comp, FFTW_MEASURE);
    c2r = fftw_plan_dft_c2r_3d(shape_x, shape_y, shape_z, field_comp, field_real,  FFTW_MEASURE);
    };
  ~RandomScalarField() {
    if (clean_switch) {
      fftw_destroy_plan(c2r);
      fftw_destroy_plan(r2c);
      fftw_free(field_real);
      fftw_free(field_comp);
      }
    #ifdef _OPENMP
          fftw_cleanup_threads();
    #else
          fftw_cleanup();
    #endif
    };

  // Fields

  int global_seed = 0;
  // for destructor
  bool clean_switch = false;
  // methods

  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  virtual double _on_grid(double* rf, fftw_complex* cf, fftw_plan forward, fftw_plan backward, const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment, const int seed) = 0;

  double on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) {
    if (recreate_fftw_plans) {
      double* field_real_temp = 0;
      fftw_complex* field_comp_temp = 0;
      field_real_temp = (double*) fftw_malloc(sizeof(double) * n[0]*n[1]*(n[2]+2));
      field_comp_temp = reinterpret_cast<fftw_complex*>(field_real_temp);
      r2c = fftw_plan_dft_r2c_3d(shape_x, shape_y, shape_z, field_real_temp, field_comp_temp, FFTW_ESTIMATE);
      c2r = fftw_plan_dft_c2r_3d(shape_x, shape_y, shape_z, field_comp_temp, field_real_temp,  FFTW_ESTIMATE);
      return _on_grid(field_real_temp, field_comp_temp, r2c, c2r, n, zeropoint, increment, seed);
      }
    else {
      return _on_grid(field_real, field_comp, r2c, c2r, n, zeropoint, increment, seed);
    }
  }


};


class RandomVectorField : public RandomField<std::vector<double>>  {
protected:
    // Fields
    int ndim = 3;
    int shape_x = 0;
    int shape_y = 0;
    int shape_z = 0;
    bool recreate_fftw_plans;
    fftw_plan c2r_x;
    fftw_plan r2c_x;
    fftw_plan c2r_y;
    fftw_plan r2c_y;
    fftw_plan c2r_z;
    fftw_plan r2c_z;
    double* field_vec_real_x = 0;
    fftw_complex* field_vec_comp_x = 0;
    double* field_vec_real_y = 0;
    fftw_complex* field_vec_comp_y = 0;
    double* field_vec_real_z = 0;
    fftw_complex* field_vec_comp_z = 0;
    // constructors
    //RandomField() : Field<G, T>() {}
    // methods

public:
  RandomVectorField() : RandomField<std::vector<double>>() {
    recreate_fftw_plans = true;
    };
  RandomVectorField(int i, int shape_x, int shape_y, int shape_z) : RandomField<std::vector<double>>(), shape_x(shape_x), shape_y(shape_y), shape_z(shape_z) {
    recreate_fftw_plans = false;
    field_vec_real_x = (double*) fftw_malloc(sizeof(fftw_complex)*shape_x*shape_y*(shape_z+2));
    field_vec_comp_x = reinterpret_cast<fftw_complex*>(field_vec_real_x);
    field_vec_real_y = (double*) fftw_malloc(sizeof(fftw_complex)*shape_x*shape_y*(shape_z+2));
    field_vec_comp_y = reinterpret_cast<fftw_complex*>(field_vec_real_y);
    field_vec_real_z = (double*) fftw_malloc(sizeof(fftw_complex)*shape_x*shape_y*(shape_z+2));
    field_vec_comp_z = reinterpret_cast<fftw_complex*>(field_vec_real_z);
    r2c_x = fftw_plan_dft_r2c_3d(shape_x, shape_y, shape_z, field_vec_real_x, field_vec_comp_x, FFTW_MEASURE);
    c2r_x = fftw_plan_dft_c2r_3d(shape_x, shape_y, shape_z, field_vec_comp_x, field_vec_real_x,  FFTW_MEASURE);
    r2c_y = fftw_plan_dft_r2c_3d(shape_x, shape_y, shape_z, field_vec_real_y, field_vec_comp_y, FFTW_MEASURE);
    c2r_y = fftw_plan_dft_c2r_3d(shape_x, shape_y, shape_z, field_vec_comp_y, field_vec_real_y,  FFTW_MEASURE);
    r2c_z = fftw_plan_dft_r2c_3d(shape_x, shape_y, shape_z, field_vec_real_z, field_vec_comp_z, FFTW_MEASURE);
    c2r_z = fftw_plan_dft_c2r_3d(shape_x, shape_y, shape_z, field_vec_comp_z, field_vec_real_z,  FFTW_MEASURE);
    };
  ~RandomVectorField() {
    if (clean_switch) {
      fftw_destroy_plan(c2r_x);
      fftw_destroy_plan(r2c_x);
      fftw_free(field_vec_real_x);
      fftw_free(field_vec_comp_x);
      fftw_destroy_plan(c2r_y);
      fftw_destroy_plan(r2c_y);
      fftw_free(field_vec_real_y);
      fftw_free(field_vec_comp_y);
      fftw_destroy_plan(c2r_z);
      fftw_destroy_plan(r2c_z);
      fftw_free(field_vec_real_z);
      fftw_free(field_vec_comp_z);
      }
    #ifdef _OPENMP
          fftw_cleanup_threads();
    #else
          fftw_cleanup();
    #endif
    };

  // for destructor
  bool clean_switch = false;
  // methods
  std::vector<double> at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
   // T c{0};
   // return c;
  };
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  double on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z) {
    throw NotImplementedException();
    //std::vector<double> c{0};
    //return c;
  };

  double on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) {
    // MISSING calculate seed!
    int nseed = 0; // needs to be updated
    if (recreate_fftw_plans) {
      field_vec_real_x = (double*) fftw_malloc(sizeof(fftw_complex)*n[0]*n[1]*(n[2]+2));
      field_vec_comp_x = reinterpret_cast<fftw_complex*>(field_vec_real_x);
      field_vec_real_y = (double*) fftw_malloc(sizeof(fftw_complex)*n[0]*n[1]*(n[2]+2));
      field_vec_comp_y = reinterpret_cast<fftw_complex*>(field_vec_real_y);
      field_vec_real_z = (double*) fftw_malloc(sizeof(fftw_complex)*n[0]*n[1]*(n[2]+2));
      field_vec_comp_z = reinterpret_cast<fftw_complex*>(field_vec_real_z);
      r2c_x = fftw_plan_dft_r2c_3d(shape_x, shape_y, shape_z, field_vec_real_x, field_vec_comp_x, FFTW_ESTIMATE);
      c2r_x = fftw_plan_dft_c2r_3d(shape_x, shape_y, shape_z, field_vec_comp_x, field_vec_real_x,  FFTW_ESTIMATE);
      r2c_y = fftw_plan_dft_r2c_3d(shape_x, shape_y, shape_z, field_vec_real_y, field_vec_comp_y, FFTW_ESTIMATE);
      c2r_y = fftw_plan_dft_c2r_3d(shape_x, shape_y, shape_z, field_vec_comp_y, field_vec_real_y,  FFTW_ESTIMATE);
      r2c_z = fftw_plan_dft_r2c_3d(shape_x, shape_y, shape_z, field_vec_real_z, field_vec_comp_z, FFTW_ESTIMATE);
      c2r_z = fftw_plan_dft_c2r_3d(shape_x, shape_y, shape_z, field_vec_comp_z, field_vec_real_z,  FFTW_ESTIMATE);
      }
    _on_grid({r2c_x, r2c_y, r2c_z}, {c2r_x, c2r_y, c2r_z}, n, zeropoint, increment, nseed);
    
  }

  void _on_grid( std::initializer_list<fftw_plan> forward, std::initializer_list<fftw_plan> backward, const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment, const int seed);
  
// this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
// original author: https://github.com/gioacchinowang
void divergence_cleaner(fftw_complex* bx, fftw_complex* by, fftw_complex* bz,  const std::vector<int> &n, const std::vector<double> &increment) const {
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < n[0]; ++i) {
      double kx = i / increment[0];
      if (i >= (n[0] + 1) / 2)
        kx -= 1. /  increment[0];
        // it's faster to calculate indices manually
      const int idx_lv1 = i * n[1] * n[2];
      for (int j = 0; j < n[1]; ++j) {
        double ky = j /  increment[1];
        if (j >= (n[1] + 1) / 2)
          ky -= 1. /  increment[1];
        const int idx_lv2 = idx_lv1 + j * n[2];
        for (int l = 0; l < n[2]/2 + 1; ++l) {
          // 0th term is fixed to zero in allocation
          if (i == 0 and j == 0 and l == 0)
            continue;
          double kz = 1. /  increment[2];
          const int idx = idx_lv2 + l;
          double k_length = 0;
          double b_length = 0;
          double b_dot_k = 0;
          std::vector<double> k{kx, ky, kz};
          b_length = static_cast<double>(bx[idx]*bx[idx] + by[idx]*by[idx] + bz[idx]*bz[idx]);
          k_length = static_cast<double>(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);
          b_dot_k = (bx[idx]*k[0] + by[idx]*k[1] + bz[idx]*k[2]);

          if (k_length == 0 or b_length == 0) {
            continue;
            }
          k_length = std::sqrt(k_length);

          const double bk_over_k = b_dot_k / k_length;
          // multiply \sqrt(3) for preserving spectral power statistically
          bx[idx] = 1.73205081 * (bx[idx] - k[0] * bk_over_k);
          by[idx] = 1.73205081 * (by[idx] - k[1] * bk_over_k);
          bz[idx] = 1.73205081 * (bz[idx] - k[2] * bk_over_k);
        } // l
      } // j
    } // i
  }
};


#endif /* RANDOMFIELD_H */
