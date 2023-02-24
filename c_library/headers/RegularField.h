#ifndef REGULARFIELD_H
#define REGULARFIELD_H


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
#include "Field.h"

class RegularScalarField : public Field<double, double*>  {
protected:
    // Fields
    const int ndim = 1;

    // Constructors

    // RegularField() : Field<T>() {};
    // Methods

public:
  // -----CONSTRUCTORS-----
  
  //using Field<double, double*> :: Field;
  ~RegularScalarField() = default;

  RegularScalarField(int shape, double zeropoint, double increment) {
    RegularScalarField : Field<double, double*>(shape=shape, zeropoint=zeropoint, increment=increment);
    class_eval = new double[array_size];
  };


  RegularScalarField(std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) {
    RegularScalarField : Field<double, double*>(grid_x, grid_y, grid_z);
    class_eval = new double[array_size];
  };
  // Fields

  // Methods

  double* on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z) {
    double* function_eval;
    evaluate_function_on_grid(function_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

  double* on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) {
    double* function_eval;
    evaluate_function_on_grid(function_eval, n, zeropoint, increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

};


class RegularVectorField : public Field<double[3], std::array<double*, 3>>  {
protected:
    // Fields

    // Constructors

    // RegularField() : Field<T>() {};
    // Methods

public:
  // Constructors
  using Field<double[3], std::array<double*, 3>> :: Field;
  virtual ~RegularVectorField() = default;

    RegularVectorField(int shape, double zeropoint, double increment) {
    RegularVectorField : Field<double, double*>(shape, zeropoint, increment);
    class_eval[0] = new double[array_size];
    class_eval[1] = new double[array_size];
    class_eval[2] = new double[array_size];
  };


  RegularVectorField(std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) {
    RegularVectorField : Field<double, double*>(grid_x, grid_y, grid_z);
    class_eval[0] = new double[array_size];
    class_eval[1] = new double[array_size];
    class_eval[2] = new double[array_size];
  };
  // Fields

  // Methods


  std::array<double*, 3> on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z) {
    std::array<double*, 3> function_eval;
    evaluate_function_on_grid(function_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

 std::array<double*, 3> on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) {
    std::array<double*, 3> function_eval;
    evaluate_function_on_grid(function_eval, n, zeropoint, increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

};

#endif
