#ifndef FIELD_H
#define FIELD_H


#include <vector>
#include <array>
#include <functional>
#include <algorithm>
#include <iostream>
#include <memory>

#include "exceptions.h"



template<typename POSTYPE, typename GRIDTYPE>
class Field {
protected:
  // -----FIELDS-----
  int ndim = 0;

  GRIDTYPE class_eval;

  bool has_grid = false;
  bool regular_grid;

  // -----CONSTRUCTORS-----

  ~Field() {};

  Field(std::array<int, 3> grid_shape, std::array<double, 3>  grid_zeropoint, std::array<double, 3>  grid_increment) : shape(grid_shape), zeropoint(grid_zeropoint), increment(grid_increment) {
    has_grid = true;
    regular_grid = true;
  };

  Field(std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : grid_x(grid_x), grid_y(grid_y), grid_z(grid_z) {
    has_grid = true;
    regular_grid = false;
    shape = {(int)grid_x.size(), (int)grid_y.size(), (int)grid_z.size()};
  };

  Field() {
    has_grid = false;
  };
  
  virtual void allocate_memory(bool not_empty, int sz) = 0;
  virtual void free_memory(bool not_empty) = 0;

public:
  // -----FIELDS-----

  std::array<int, 3> shape;
  std::array<double, 3> zeropoint;
  std::array<double, 3> increment;
  std::vector<double> grid_x;
  std::vector<double> grid_y;
  std::vector<double> grid_z;

  // -----METHODS-----

  size_t array_size() {
    if (has_grid) {
      size_t arsz = 1;
      for (const int sh : shape) 
        arsz = arsz*sh;
      return arsz;
      }
    else 
      return 0;
  }
  // -----Interface functions-----

  // Evaluate the model at Galactic position (x, y, z). 
  virtual POSTYPE at_position(const double &x, const double &y, const double &z) const = 0;

   // Evaluate the model on a grid. 
   //The grid may be provided as three vectors containting the x, y and z coordinates (useful for irregular grids) or via providing the zeropoint, increment and number of pixels along each axis (thereby defining a regular grid). 
   // Both options maybe provided in the initialization of the class (thereby potentially ) 
  
  virtual GRIDTYPE on_grid(const int seed = 0) = 0;

  virtual GRIDTYPE on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, const int seed = 0) = 0;

  virtual GRIDTYPE on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed = 0) = 0;

  // This is the interface function to CRPRopa
  POSTYPE getField(const std::array<double, 3> &pos_arr) const {
    return at_position(pos_arr[0], pos_arr[1], pos_arr[2]);
  }
  
  // -----Helper functions-----

// Evaluate scalar valued functions on irregular grids
  void evaluate_function_on_grid(double *fval, const std::vector<double> &ggx, const std::vector<double> &ggy, const std::vector<double> &ggz, std::function<double(double, double, double)> eval) {
    std::vector<int> size{(int)ggx.size(), (int)ggy.size(), (int)ggz.size()};
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
  void evaluate_function_on_grid(double *fval, const std::array<int, 3> &size, const std::array<double, 3> &zp, const std::array<double, 3> &inc, std::function<double(double, double, double)> eval) {
     for (int i=0; i < size[0]; i++) {
         int m = i*size[0];
         for (int j=0; j < size[1]; j++) {
             int n = j*size[1];
             for (int k=0; k < size[2]; k++) {
                 fval[m + n + k] = eval(zp[0] + i*inc[0], zp[1] + j*inc[1], zp[2] + k*inc[2]);
         }   }   }
    }


  // Evaluate vector valued functions on irregular grids (replace with template for ndim?)
  void evaluate_function_on_grid(std::array<double*, 3>  fval, const std::vector<double> &ggx, const std::vector<double> &ggy, const std::vector<double> &ggz, std::function<std::array<double, 3>(double, double, double)> eval) {
    std::vector<int> size{(int)ggx.size(), (int)ggy.size(), (int)ggz.size()};
    std::cout << "The ref to b_grid when filled " << &fval << " \n\n";
     for (int i=0; i < size[0]; i++) {
         int m = i*size[1]*size[2];
         for (int j=0; j < size[1]; j++) {
             int n = j*size[2];
             for (int k=0; k < size[2]; k++) {
                 std::array<double, 3> v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
                 fval[0][m + n + k] = v[0];
                 fval[1][m + n + k] = v[1];
                 fval[2][m + n + k] = v[2];
              }
          }   
      }   
   }

  // Evaluate vector valued functions on regular grids
  void evaluate_function_on_grid(std::array<double*, 3>  fval, const std::array<int, 3> &size,
                                  const std::array<double, 3> &zp, const std::array<double, 3> &inc,
                                  std::function<std::array<double, 3>(double, double, double)> eval) {
      for (int i=0; i < size[0]; i++) {
          int m = i*size[1]*size[2];
          for (int j=0; j < size[1]; j++) {
              int n = j*size[2];
              for (int k=0; k < size[2]; k++) {
                  std::array<double, 3> v = eval(zp[0] + i*inc[0], zp[1] + j*inc[1], zp[2] + k*inc[2]);
                  fval[0][m + n + k] = v[0];
                  fval[1][m + n + k] = v[1];
                  fval[2][m + n + k] = v[2];
                  
          }   }   }
    }

};

#endif