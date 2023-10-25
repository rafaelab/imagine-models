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
  bool has_fftw_wisdom = false;

public:
  // Constructors
  using Field<POSTYPE, GRIDTYPE> :: Field;

  // Fields

  bool apply_spectrum = true;


  // methods
  POSTYPE at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
    // Here comes the interpolator
  }
  
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  virtual double calculate_fourier_sigma(const double &abs_k, const double &dk) const = 0;

    // This function is the place where the global routine should be implemented, i.e. how the spatial profile modifies the random field, and if divergence cleaning needs to be performed. 
  virtual void _on_grid(GRIDTYPE val, const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) = 0;

  void draw_random_numbers_complex(fftw_complex* vec,  const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed); 

  GRIDTYPE on_grid(const std::vector<double>  &grid_x, const std::vector<double>  &grid_y, const std::vector<double>  &grid_z, const int seed) {
    throw NotImplementedException();
  }

  void remove_padding(double* val, const std::array<int, 3> &shp, const int pad);

  double simple_spectrum(const double &abs_k, const double &dk, const double &k0, const double &s) const {
      double pi = 3.141592653589793;
      const double unit = 1. / (4 * pi * abs_k * abs_k); 
      const double dP = unit / std::pow(abs_k + k0, s);
      //double norm = 1. / ((s - 1.) * std::pow(k0, (s - 1))); // normalize to unity
      return dP; // / norm;
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


class RandomScalarField : public RandomField<number, double*>  {
protected:
    // Fields
    fftw_plan c2r;
    fftw_plan r2c;
    // fftw_complex* grid_eval_comp = NULL;
    // constructors
    //RandomField() : Field<G, T>() {}
    // methods

  double* allocate_memory(std::array<int, 3> shp) {
    int newshp2;
    if (shp[2] % 2) {
      newshp2 = shp[2] + 1;
    }
    else {
      newshp2 = shp[2] + 2;
    }
    double* grid_eval = (double*) fftw_alloc_real(shp[0]*shp[1]*newshp2);
    return grid_eval;  
  }  

  void free_memory(double* grid_eval) {
      fftw_free(grid_eval);
  }

  fftw_complex* construct_plans(double* grid_eval, std::array<int, 3> shp) {
    fftw_complex* grid_eval_comp = reinterpret_cast<fftw_complex*>(grid_eval);
    if (has_fftw_wisdom) {
        const char *filename = "ImagineModelsRandomFieldScalar"; // FIX this, currently overwrites
        int fftw_im_wisdom_to_filename(*filename);
    }
    r2c = fftw_plan_dft_r2c_3d(shp[0], shp[1], shp[2], grid_eval, grid_eval_comp, FFTW_ESTIMATE);
    c2r = fftw_plan_dft_c2r_3d(shp[0], shp[1], shp[2], grid_eval_comp, grid_eval,  FFTW_ESTIMATE);
    created_fftw_plans = true;
    return grid_eval_comp;
  }

  void destroy_plans() {
    if (created_fftw_plans) {
      fftw_destroy_plan(c2r);
      fftw_destroy_plan(r2c);
    }

  }

public:
  RandomScalarField() : RandomField<number, double*>() {
    };

  RandomScalarField(std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  increment) : RandomField<number, double*>(shape, reference_point, increment) {
    //accumulate wisdom
    double* grid_eval = allocate_memory(shape);
    fftw_complex* grid_eval_comp = reinterpret_cast<fftw_complex*>(grid_eval);
    fftw_plan r2c_temp = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], grid_eval, grid_eval_comp, FFTW_MEASURE);
    fftw_plan c2r_temp = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], grid_eval_comp, grid_eval, FFTW_MEASURE);
    //save wisdom
    const char *filename = "ImagineModelsRandomScalarField";
    int fftw_export_wisdom_to_filename(*filename);
    //destroy wisdom
    has_fftw_wisdom = true;
    fftw_destroy_plan(c2r_temp);
    fftw_destroy_plan(r2c_temp); //delete plans, but wisdom will be kept!
    fftw_free(grid_eval);
    };

  ~RandomScalarField() {
    //delete wisdom
    destroy_plans();
    fftw_forget_wisdom();
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


  double* on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) {
    if (initialized_with_grid) 
      throw GridException();
    double* grid_eval = allocate_memory(shp);
    _on_grid(grid_eval, shp, rpt, inc, seed);
    return grid_eval;
  }

  double* on_grid(const int seed) {
    if (not initialized_with_grid) 
      throw GridException();
    double* grid_eval = allocate_memory(internal_shape);
    const char *filename = "ImagineModelsRandomScalarField";
    int fftw_import_wisdom_from_filename(*filename);
    _on_grid(grid_eval, internal_shape, internal_ref_point, internal_increment, seed);
    return grid_eval;
  }
};



class RandomVectorField : public RandomField<vector, std::array<double*, 3>>  {
protected:
    // Fields

  std::array<fftw_plan, 3> r2c;
  std::array<fftw_plan, 3> c2r;

  std::array<double*, 3> allocate_memory(std::array<int, 3> shp) {
    std::array<double*, 3> grid_eval;
    int newshp2;
    if (shp[2] % 2) {
      newshp2 = shp[2] + 1;
    }
    else {
      newshp2 = shp[2] + 2;
    }
    for (int i=0; i < ndim; ++i) {
      grid_eval[i] = (double*) fftw_alloc_real(shp[0]*shp[1]*newshp2);
      }
    return grid_eval;  
  }  

  void free_memory(std::array<double*, 3> grid_eval) {
    for (int i=0; i < ndim; ++i) { 
      fftw_free(grid_eval[i]);
    }
  }

  std::array<fftw_complex*, 3> construct_plans(std::array<double*, 3> grid_eval, std::array<int, 3> shp) {

  std::array<fftw_complex*, 3> grid_eval_comp;
    for (int i=0; i < ndim; ++i) {
      grid_eval_comp[i] = reinterpret_cast<fftw_complex*>(grid_eval[i]);
      r2c[i] = fftw_plan_dft_r2c_3d(shp[0], shp[1], shp[2], grid_eval[i], grid_eval_comp[i], FFTW_ESTIMATE);
      c2r[i] = fftw_plan_dft_c2r_3d(shp[0], shp[1], shp[2], grid_eval_comp[i], grid_eval[i],  FFTW_ESTIMATE);
    }
    created_fftw_plans = true;
  return grid_eval_comp;
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

  RandomVectorField(std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  increment) : RandomField(shape, reference_point, increment) {
    int newshp2;
    if (shape[2] % 2) {
      newshp2 = shape[2] + 1;
    }
    else {
      newshp2 = shape[2] + 2;
    }
    double* val_temp = (double*) fftw_alloc_real(shape[0]*shape[1]*newshp2);
    fftw_complex* val_temp_comp = reinterpret_cast<fftw_complex*>(val_temp);
    fftw_plan r2c_temp = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], val_temp, val_temp_comp, FFTW_MEASURE);
    fftw_plan c2r_temp = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], val_temp_comp, val_temp,  FFTW_MEASURE);
    //save wisdom
    const char *filename = "ImagineModelsRandomVectorField";
    int fftw_export_wisdom_to_filename(*filename);
    has_fftw_wisdom = true;
    fftw_free(val_temp);
    fftw_destroy_plan(c2r_temp);
    fftw_destroy_plan(r2c_temp); // plans are destroyed but wisdom is saved!
    };

  ~RandomVectorField() {
    //delete wisdom
    destroy_plans();
    fftw_forget_wisdom();
    #ifdef _OPENMP
          fftw_cleanup_threads();
    #else
          fftw_cleanup();
    #endif
    };

  // fields
  const int ndim = 3;

  bool clean_divergence = true;  
  // methods
  vector at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
   // T c{0};
   // return c;
  }
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;


  std::array<double*, 3> on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) {
    if (initialized_with_grid) 
      throw GridException();
    std::array<double*, 3> grid_eval = allocate_memory(shp);
    _on_grid(grid_eval, shp, rpt, inc, seed);
    return grid_eval;
    
  }

  std::array<double*, 3> on_grid(const int seed) {
    if (not initialized_with_grid) 
      throw GridException();
    std::array<double*, 3> grid_eval = allocate_memory(internal_shape);
    const char *filename = "ImagineModelsRandomVectorField";
    int fftw_import_wisdom_from_filename(*filename);
    _on_grid(grid_eval, internal_shape, internal_ref_point, internal_increment, seed);
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
