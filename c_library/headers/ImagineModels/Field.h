#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <array>
#include <functional>
#include <algorithm>
#include <iostream>
#include <memory>

#include "exceptions.h"

#if autodiff_FOUND
    #include <autodiff/forward/real.hpp>
    #include <autodiff/forward/dual.hpp>
    #include <autodiff/forward/real/eigen.hpp>
    namespace ad = autodiff;
    typedef ad::real number;
    typedef ad::VectorXreal vector;
#else
    typedef double number;  // only used for differentiable numbers! 
    typedef std::array<double, 3> vector;
#endif


template<typename POSTYPE, typename GRIDTYPE>
class Field {
protected:
  // -----FIELDS-----
  int ndim = 0;

  bool initialized_with_grid;
  bool regular_grid;
  

  // -----CONSTRUCTORS-----

  Field(std::array<int, 3> shape, std::array<double, 3>  ref_point, std::array<double, 3>  increment) : internal_shape(shape), internal_ref_point(ref_point), internal_increment(increment) {
    initialized_with_grid = true;
    regular_grid = true;
  };

  Field(std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : internal_grid_x(grid_x), internal_grid_y(grid_y), internal_grid_z(grid_z) {
    initialized_with_grid = true;
    regular_grid = false;
    internal_shape = {(int)grid_x.size(), (int)grid_y.size(), (int)grid_z.size()};
  };

  Field() {
    initialized_with_grid = false;
  };
  
  virtual GRIDTYPE allocate_memory(std::array<int, 3> shp) = 0;
  virtual void free_memory(GRIDTYPE grid_eval) = 0;

public:

  ~Field() {};
  // -----FIELDS-----

  std::array<int, 3> internal_shape;
  std::array<double, 3> internal_ref_point;
  std::array<double, 3> internal_increment;
  std::vector<double> internal_grid_x;
  std::vector<double> internal_grid_y;
  std::vector<double> internal_grid_z;

  // -----METHODS-----

  // -----Interface functions-----

  // Evaluate the model at Galactic position (x, y, z). 
  virtual POSTYPE at_position(const double &x, const double &y, const double &z) const = 0;

   // Evaluate the model on a grid. 
   //The grid may be provided as three vectors containting the x, y and z coordinates (useful for irregular grids) or via providing the zeropoint, increment and number of pixels along each axis (thereby defining a regular grid). 
   // Each of the options maybe provided in the initialization of the class (thereby potentially allowing precomputation) 
  
  virtual GRIDTYPE on_grid(const int seed = 0) = 0;

  virtual GRIDTYPE on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, const int seed = 0) = 0;

  virtual GRIDTYPE on_grid(const std::array<int, 3> &shape, const std::array<double, 3> &ref_point, const std::array<double, 3> &increment, const int seed = 0) = 0;
  
  // -----Helper functions-----

  // Get number of pixels
  size_t grid_size(std::array<int, 3> shp) {
      size_t gsz = 1;
      for (const int &sh : shp) {
        gsz = gsz*sh;}
      return gsz;
  }


  // Evaluate scalar valued functions on regular grids
  template <typename FRTYPE> 
  void evaluate_function_on_grid(double *fval, const std::array<int, 3> &size, const std::array<double, 3> &rpt,  const std::array<double, 3> &inc, std::function<FRTYPE(double, double, double)> eval) {
    for (int i=0; i < size[0]; i++) {
      int m = i*size[0];
      for (int j=0; j < size[1]; j++) {
        int n = j*size[1];
        for (int k=0; k < size[2]; k++) {
          FRTYPE v = eval(rpt[0] + i*inc[0], rpt[1] + j*inc[1], rpt[2] + k*inc[2]);
                    fval[m + n + k] = static_cast<double>(v);
        }   
      }   
    }
  } 

// Evaluate scalar valued functions on irregular grids
  template <typename FRTYPE> 
  void evaluate_function_on_grid(double *fval, const std::vector<double> &ggx, const std::vector<double> &ggy, const std::vector<double> &ggz, std::function<FRTYPE(double, double, double)> eval) {
    std::vector<int> size{(int)ggx.size(), (int)ggy.size(), (int)ggz.size()};
    for (int i=0; i < size[0]; i++) {
      int m = i*size[1]*size[2];
      for (int j=0; j < size[1]; j++) {
        int n = j*size[2];
        for (int k=0; k < size[2]; k++) {
          FRTYPE v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
          fval[m + n + k] = static_cast<double>(v);
        }   
      }   
    }
  }


  // Evaluate vector valued functions on irregular grids (replace with template for ndim?)
  template <typename FRTYPE> 
  void evaluate_function_on_grid(std::array<double*, 3>  fval, const std::vector<double> &ggx, const std::vector<double> &ggy, const std::vector<double> &ggz, std::function<FRTYPE(double, double, double)> eval) {
    std::vector<int> size{(int)ggx.size(), (int)ggy.size(), (int)ggz.size()};
    for (int i=0; i < size[0]; i++) {
      int m = i*size[1]*size[2];
      for (int j=0; j < size[1]; j++) {
        int n = j*size[2];
        for (int k=0; k < size[2]; k++) {
          FRTYPE v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
          fval[0][m + n + k] = static_cast<double>(v[0]);
          fval[1][m + n + k] = static_cast<double>(v[1]);
          fval[2][m + n + k] = static_cast<double>(v[2]);
        }      
      }   
    }   
  }

  // Evaluate vector valued functions on regular grids
  template <typename FRTYPE> 
  void evaluate_function_on_grid(std::array<double*, 3>  fval, const std::array<int, 3> &size, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, std::function<FRTYPE(double, double, double)> eval) {
    for (int i=0; i < size[0]; i++) {
      int m = i*size[1]*size[2];
      for (int j=0; j < size[1]; j++) {
        int n = j*size[2];
        for (int k=0; k < size[2]; k++) {
          FRTYPE v = eval(rpt[0] + i*inc[0], rpt[1] + j*inc[1], rpt[2] + k*inc[2]);
          fval[0][m + n + k] = static_cast<double>(v[0]);
          fval[1][m + n + k] = static_cast<double>(v[1]);
          fval[2][m + n + k] = static_cast<double>(v[2]);        
        }   
      }  
    }
  }

};

#endif