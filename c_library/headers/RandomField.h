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

public:
  // Constructors
  using Field<POSTYPE, GRIDTYPE> :: Field;

  // Fields
  int seed = 0;
  double rms;
  double k0;
  double k1;
  double a0;
  double a1;

  // methods
  POSTYPE at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
    // Here comes the interpolator
    POSTYPE c{0};
    return c;
  }
  
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  GRIDTYPE on_grid(const std::vector<double>  &grid_x, const std::vector<double>  &grid_y, const std::vector<double>  &grid_z) {
    throw NotImplementedException();
    GRIDTYPE c{0};
    return c;
  }

  // below breaks for the scalar field
  void draw_random_numbers(std::array<fftw_complex*, 3> c_fields, const std::array<int, 3> &grid_shape, const std::array<double, 3>  &grid_increment, const int seed);

  double spectrum(const double &abs_k, const double &rms, const double &k0, const double &k1, const double &a0, const double &a1) const;

};

//include the random field method implementations (at this position, due to the use of templates)
#include "RandomField.tpp"

class RandomScalarField : public RandomField<double, double*>  {
protected:
    // Fields
    bool recreate_fftw_plans;
    fftw_plan c2r;
    fftw_plan r2c;
    fftw_complex* class_eval_comp = 0;
    // constructors
    //RandomField() : Field<G, T>() {}
    // methods

  void allocate_memory(bool not_empty, int arr_sz) {
    if (not_empty) {
    class_eval = (double*) fftw_malloc(sizeof(fftw_complex) * shape[0]*shape[1]*(shape[2]+2));
     class_eval_comp = reinterpret_cast<fftw_complex*>( class_eval);
    r2c = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], class_eval, class_eval_comp, FFTW_MEASURE);
    c2r = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], class_eval_comp, class_eval,  FFTW_MEASURE);
      }
    else {
      class_eval = 0;
      };
    }

  void free_memory(bool not_empty) {
    if (not_empty) {
      fftw_free(class_eval);
      fftw_destroy_plan(c2r);
      fftw_destroy_plan(r2c);
    };
  }

public:
  RandomScalarField() : RandomField<double, double*>() {
    };

  RandomScalarField(std::array<int, 3>  shape, std::array<double, 3>  zeropoint, std::array<double, 3>  increment) : RandomField<double, double*>(shape, zeropoint, increment) {
    allocate_memory(has_grid, array_size);

    };

  ~RandomScalarField() {
    free_memory(has_grid);
    #ifdef _OPENMP
          fftw_cleanup_threads();
    #else
          fftw_cleanup();
    #endif
    };

  // Fields

  const int ndim = 1;
  // for destructor
  bool clean_switch = false;
  // methods

  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  virtual double* _on_grid(double* rf, fftw_complex* cf, fftw_plan forward, fftw_plan backward, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) = 0;

  double* on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {
    double* field_real_temp = 0;
    fftw_complex* field_comp_temp = 0;
    field_real_temp = (double*) fftw_malloc(sizeof(fftw_complex) * shape[0]*shape[1]*(shape[2]+2));
    field_comp_temp = reinterpret_cast<fftw_complex*>( field_real_temp);
    fftw_plan r2c_temp = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], field_real_temp,field_comp_temp, FFTW_MEASURE);
    fftw_plan c2r_temp = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], field_comp_temp, field_real_temp,  FFTW_MEASURE);
    return _on_grid(field_real_temp, field_comp_temp, r2c_temp, c2r_temp, grid_shape, grid_zeropoint, grid_increment, seed);
  }

  double* on_grid(const int seed) {
    return _on_grid(class_eval, class_eval_comp, r2c, c2r, shape, zeropoint, increment, seed);
  }


};


class RandomVectorField : public RandomField<std::array<double, 3>, std::array<double*, 3>>  {
protected:
    // Fields

    std::array<fftw_plan, 3> r2c;
    std::array<fftw_plan, 3> c2r;
    std::array<fftw_complex*, 3> class_eval_comp;

    // constructors
    //RandomField() : Field<G, T>() {}
    // methods


  void allocate_memory(bool not_empty, int arr_sz) {
    if (not_empty) {
      for (int i=0; i < ndim; ++i) {
        class_eval[i] = (double*) fftw_malloc(sizeof(fftw_complex) * shape[0]*shape[1]*(shape[2]+2));
        class_eval_comp[i]= reinterpret_cast<fftw_complex*>( class_eval[i]);
        r2c[i] = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], class_eval[i], class_eval_comp[i], FFTW_MEASURE);
        c2r[i] = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], class_eval_comp[i], class_eval[i],  FFTW_MEASURE);
        } 
      }
    else {
      class_eval = {0, 0, 0};
      };
    }

  void free_memory(bool not_empty) {
    if (not_empty) {
      for (int i=0; i < ndim; ++i) { 
        fftw_free(class_eval[i]);
        fftw_destroy_plan(c2r[i]);
        fftw_destroy_plan(r2c[i]);
      }
    };
  }
public:


  RandomVectorField() : RandomField() {
    allocate_memory(has_grid, array_size);
    };
  RandomVectorField(std::array<int, 3>  shape, std::array<double, 3>  zeropoint, std::array<double, 3>  increment) : RandomField(shape, zeropoint, increment) {
    allocate_memory(has_grid, array_size);
    };
  ~RandomVectorField() {
    free_memory(has_grid);
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

  std::array<double*, 3> on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, int seed = 0) {
    throw NotImplementedException();
    //std::vector<double> c{0};
    //return c;
  }

  std::array<double*, 3> on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, int nseed) {
    std::array<fftw_plan, 3> r2c_temp;
    std::array<fftw_plan, 3> c2r_temp;
    std::array<fftw_complex*, 3> eval_comp_temp;
    std::array<double*, 3> eval_temp;
    for (int i=0; i < ndim; ++i) {
      eval_temp[i] = (double*) fftw_malloc(sizeof(fftw_complex) * shape[0]*shape[1]*(shape[2]+2));
      eval_comp_temp[i]= reinterpret_cast<fftw_complex*>( eval_temp[i]);
      r2c[i] = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], eval_temp[i], eval_comp_temp[i], FFTW_MEASURE);
      c2r[i] = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], eval_comp_temp[i], eval_temp[i],  FFTW_MEASURE);

    } 
    _on_grid(eval_temp, eval_comp_temp, r2c_temp, c2r_temp, grid_shape, grid_zeropoint, grid_increment, nseed);
    return eval_temp;
    
  }

  std::array<double*, 3> on_grid(int nseed) {
    _on_grid(class_eval, class_eval_comp, r2c, c2r, shape, zeropoint, increment, nseed);
    return class_eval;
  }
  

  virtual void _on_grid(std::array<double*, 3> rf, std::array<fftw_complex*, 3> cf, std::array<fftw_plan, 3> forward, std::array<fftw_plan, 3> backward, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) = 0;
  
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
