#ifndef FIELD_H
#define FIELD_H


#include <vector>
#include <cassert>
#include <stdexcept>
#include <functional>
#include <iostream>
#include <algorithm>
#include <fftw3.h>
#include <random>
#include <memory>

#include "exceptions.h"



template<typename T>
class Field {
protected:
    // Fields

    // Constructors
    Field() {};
    //Field(const Field<G, T>& f) {};
    // Methods

public:
  // Constructors
  virtual ~Field() = default;
  // Fields

  // Methods
  virtual T at_position(const double &x, const double &y, const double &z) const = 0;
  virtual double on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z) = 0;
  virtual double on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) = 0;

  // This is the interface function to CRPRopa
  T getField(const std::vector<double> &pos_vec) const {
    return at_position(pos_vec[0], pos_vec[1], pos_vec[2]);
  }
    // Evaluate scalar valued functions on irregular grids
  void evaluate_function_on_grid(double *b, const std::vector<int> &size,
                                 const double &ggx, const double &ggy, const double &ggz,
                                 std::function<double(double, double, double)> eval) {
     for (int i=0; i < size[0]; i++) {
         int m = i*size[1]*size[2];
         for (int j=0; j < size[1]; j++) {
             int n = j*size[2];
             for (int k=0; k < size[2]; k++) {
                 b[m + n + k] = eval(ggx[i], ggy[j], ggz[k]);
         }   }   }
     }


    // Evaluate vector valued functions on irregular grids
  void evaluate_function_on_grid(double* bx, double* by, double* bz, const std::vector<int> &size,
                                 const std::vector<double> &ggx, const std::vector<double> &ggy, const std::vector<double> &ggz,
                                 std::function<std::vector<double>(double, double, double)> eval) {
     for (int i=0; i < size[0]; i++) {
         int m = i*size[1]*size[2];
         for (int j=0; j < size[1]; j++) {
             int n = j*size[2];
             for (int k=0; k < size[2]; k++) {
                 std::vector<double> v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
                 bx[m + n + k] = v[0];
                 by[m + n + k] = v[1];
                 bz[m + n + k] = v[2];
              }
          }   
      }   
   }

     // Evaluate vector valued functions on regular grids
   void evaluate_function_on_grid(double* bx, double* by, double* bz, const std::vector<int> &size,
                                  const std::vector<double> &zp, const std::vector<double> &inc,
                                  std::function<std::vector<double>(double, double, double)> eval) {
      for (int i=0; i < size[0]; i++) {
          int m = i*size[1]*size[2];
          for (int j=0; j < size[1]; j++) {
              int n = j*size[2];
              for (int k=0; k < size[2]; k++) {
                  std::vector<double> v = eval(zp[0] + i*inc[0], zp[1] + j*inc[1], zp[2] + k*inc[2]);
                  bx[m + n + k] = v[0];
                  by[m + n + k] = v[1];
                  bz[m + n + k] = v[2];
                  
          }   }   }
    }

    // Evaluate scalar valued functions on regular grids
  void evaluate_function_on_grid(double *b, const std::vector<int> &size,
                                 const std::vector<double> &zp, const std::vector<double> &inc,
                                 std::function<double(double, double, double)> eval) {
     for (int i=0; i < size[0]; i++) {
         int m = i*size[0];
         for (int j=0; j < size[1]; j++) {
             int n = j*size[1];
             for (int k=0; k < size[2]; k++) {
                 b[m + n + k] = eval(zp[0] + i*inc[0], zp[1] + j*inc[1], zp[2] + k*inc[2]);
         }   }   }
     }

};

#endif