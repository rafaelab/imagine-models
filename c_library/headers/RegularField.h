#ifndef REGULARFIELD_H
#define REGULARFIELD_H


#include <vector>
#include <array>
#include <stdexcept>
#include <functional>
#include <iostream>
#include <algorithm>

#include "exceptions.h"
#include "Field.h"

class RegularScalarField : public Field<double, double*>  {
protected:
    // Fields

    // Constructors

    // RegularField() : Field<T>() {};
    // Methods
void allocate_memory(bool not_empty, int arr_sz) {
  if (not_empty) {
    class_eval = new double[arr_sz];
    }
  else {
    class_eval = 0;
    };
  }

void free_memory(bool not_empty) {
  if (not_empty) {
    delete class_eval;
  };
}

public:
  // -----CONSTRUCTORS-----
  
  using Field<double, double*> :: Field;

  // Fields
  const int ndim = 1;
  // Methods

  double* on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z) {
    double* function_eval;
    evaluate_function_on_grid(function_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

  double* on_grid(const std::array<int, 3> &n, const std::array<double, 3> &zeropoint, const std::array<double, 3> &increment) {
    double* function_eval;
    evaluate_function_on_grid(function_eval, n, zeropoint, increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

};


class RegularVectorField : public Field<std::array<double, 3>, std::array<double*, 3>>  {
protected:
    // Fields

    // Constructors

    // RegularField() : Field<T>() {};
    // Methods
void allocate_memory(bool not_empty, int arr_sz) {
  if (not_empty) {
    class_eval[0] = new double[arr_sz];
    class_eval[1] = new double[arr_sz];
    class_eval[2] = new double[arr_sz];
    }
  else {
    class_eval[0] = 0;
    class_eval[1] = 0;
    class_eval[2] = 0;
    };
  }

void free_memory(bool not_empty) {
  if (not_empty) {
    delete class_eval[0];
    delete class_eval[1];
    delete class_eval[2];
  };
}


public:

  // Constructors
  using Field<std::array<double, 3>, std::array<double*, 3>> :: Field;
  // Fields

  // Methods


  std::array<double*, 3> on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z) {
    std::array<double*, 3> function_eval;
    evaluate_function_on_grid(function_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

 std::array<double*, 3> on_grid(const std::array<int, 3> &n, const std::array<double, 3> &zeropoint, const std::array<double, 3> &increment) {
    std::array<double*, 3> function_eval;
    evaluate_function_on_grid(function_eval, n, zeropoint, increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

};

#endif
