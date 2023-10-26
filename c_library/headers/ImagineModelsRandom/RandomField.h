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



#endif /* RANDOMFIELD_H */
