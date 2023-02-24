#ifndef FIELD_H
#define FIELD_H


#include <vector>
#include <array>
#include <functional>
#include <algorithm>
#include <memory>

#include "exceptions.h"



template<typename POSTYPE, typename GRIDTYPE>
class Field {
protected:
  // -----FIELDS-----
  int ndim;

  int[3] shape;
  double[3] zeropoint;
  double[3] increment;
  int array_size;

  static std::vector<double> grid_x;
  static std::vector<double> grid_y;
  static std::vector<double> grid_z;

  GRIDTYPE class_eval;

  bool has_grid;
  bool regular_grid;

  // -----CONSTRUCTORS-----

  ~Field() = default;

  Field(int shape, double zeropoint, double increment) : shape(shape), zeropoint(zeropoint), increment(increment) {
    array_size = 1;
    for (const int sh : shape) 
      array_size = array_size*sh;
    has_grid = true;
    regular_grid = true;
  };

  Field(std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : grid_x(grid_x), grid_y(grid_y), grid_z(grid_z) {
    has_grid = true;
    regular_grid = false;
    shape = {(int)grid_x.size(), (int)grid_y.size(), (int)grid_z.size()};
    array_size = 1;
    for (const int sh : shape) 
      array_size = array_size*sh;
  };

  Field() {
    has_grid = false;
  };


public:

  // -----METHODS-----

  // -----Interface functions-----

  // Evaluate the model at Galactic position (x, y, z). 
  virtual POSTYPE at_position(const double &x, const double &y, const double &z) const = 0;

   // Evaluate the model on a grid. 
   //The grid may be provided as three vectors containting the x, y and z coordinates (useful for irregular grids) or via providing the zeropoint, increment and number of pixels along each axis (thereby defining a regular grid). 
   // Both options maybe provided in the initialization of the class (thereby potentially ) 
  
  GRIDTYPE on_grid() {
    if (has_grid) { 
      if (regular_grid) {
        evaluate_function_on_grid(class_eval, shape, zeropoint, increment, at_position);
        return class_eval;
      }
      else {
        evaluate_function_on_grid(class_eval, grid_x, grid_y, grid_z, at_position);
        return class_eval;
      }
    }
    else
      throw GridException();
  }

  virtual GRIDTYPE on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z) = 0;

  virtual GRIDTYPE on_grid(const std::array<int> &n, const std::array<double> &zeropoint, const std::array<double> &increment) = 0;

  // This is the interface function to CRPRopa
  double getField(const std::vector<double> &pos_vec) const {
    return at_position(pos_vec[0], pos_vec[1], pos_vec[2]);
  }
  
  // -----Helper functions-----

// Evaluate scalar valued functions on irregular grids
  static void evaluate_function_on_grid(double *fval, const std::vector<double> &ggx, const std::vector<double> &ggy, const std::vector<double> &ggz, std::function<double(double, double, double)> eval) {
    std::vector<int> size{(int)ggx.size(), ggy.size(), ggz.size()};
    for (int i=0; i < size[0]; i++) {
      int m = i*size[1]*size[2];
      for (int j=0; j < size[1]; j++) {
        int n = j*size[2];
        for (int k=0; k < size[2]; k++) {
          fval[m + n + k] = eval(ggx[i], ggy[j], ggz[k]);
          }   
        }   
      }
    }

  // Evaluate scalar valued functions on regular grids
  static void evaluate_function_on_grid(double *fval, const std::vector<int> &size, const std::vector<double> &zp, const std::vector<double> &inc, std::function<double(double, double, double)> eval) {
     for (int i=0; i < size[0]; i++) {
         int m = i*size[0];
         for (int j=0; j < size[1]; j++) {
             int n = j*size[1];
             for (int k=0; k < size[2]; k++) {
                 fval[m + n + k] = eval(zp[0] + i*inc[0], zp[1] + j*inc[1], zp[2] + k*inc[2]);
         }   }   }
     }


  // Evaluate vector valued functions on irregular grids (replace with template for ndim?)
  static void evaluate_function_on_grid(std::array<double*, 3>  fval, const std::vector<double> &ggx, const std::vector<double> &ggy, const std::vector<double> &ggz, std::function<double[3](double, double, double)> eval) {
    std::vector<int> size{(int)ggx.size(), ggy.size(), ggz.size()};
     for (int i=0; i < size[0]; i++) {
         int m = i*size[1]*size[2];
         for (int j=0; j < size[1]; j++) {
             int n = j*size[2];
             for (int k=0; k < size[2]; k++) {
                 T v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
                 fval[0][m + n + k] = v[0];
                 fval[1][m + n + k] = v[1];
                 fval[2][m + n + k] = v[2];
              }
          }   
      }   
   }

  // Evaluate vector valued functions on regular grids
  static void evaluate_function_on_grid(std::array<double*, 3>  fval, const std::vector<int> &size,
                                  const std::vector<double> &zp, const std::vector<double> &inc,
                                  std::function<std::vector<double>(double, double, double)> eval) {
      for (int i=0; i < size[0]; i++) {
          int m = i*size[1]*size[2];
          for (int j=0; j < size[1]; j++) {
              int n = j*size[2];
              for (int k=0; k < size[2]; k++) {
                  std::vector<double> v = eval(zp[0] + i*inc[0], zp[1] + j*inc[1], zp[2] + k*inc[2]);
                  fval[0][m + n + k] = v[0];
                  fval[1][m + n + k] = v[1];
                  fval[2][m + n + k] = v[2];
                  
          }   }   }
    }

};

#endif