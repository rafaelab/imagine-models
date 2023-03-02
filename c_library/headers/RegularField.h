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

  ~RegularScalarField() {free_memory(has_grid);};
  
  RegularScalarField () : Field<double, double*>() {
      size_t ar_sz = array_size();
      allocate_memory(has_grid, ar_sz);
      };

  RegularScalarField (std::array<int, 3>  shape, std::array<double, 3>  zeropoint, std::array<double, 3>  increment) : Field<double, double*>(shape, zeropoint, increment) {
      size_t ar_sz = array_size();
      allocate_memory(has_grid, ar_sz);
      };

  RegularScalarField (std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : Field<double, double*>(grid_x, grid_y, grid_z) {
      size_t ar_sz = array_size();
      allocate_memory(has_grid, ar_sz);
      };

  // Fields
  const int ndim = 1;
  // Methods

  double* on_grid(int seed = 0) {
    if (not has_grid) {
      throw GridException();
    }
    if (has_grid) { 
      if (regular_grid) {
        evaluate_function_on_grid(class_eval, shape, zeropoint, increment, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
        return class_eval;
      }
      else {
        evaluate_function_on_grid(class_eval, grid_x, grid_y, grid_z, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
        return class_eval;
      }
    }
    else
      return class_eval;

  }

  double* on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, const int seed = 0) {
    double *function_eval = new double[grid_x.size()*grid_y.size()*grid_z.size()];
    //auto function_eval = std::make_shared<double>();
    evaluate_function_on_grid(function_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

  double* on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed = 0) {
    double *function_eval = new double[grid_x.size()*grid_y.size()*grid_z.size()];
    evaluate_function_on_grid(function_eval, grid_shape, grid_zeropoint, grid_increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

};


class RegularVectorField : public Field<std::array<double, 3>, std::array<double*, 3>>  {
protected:
    // Fields

    // RegularField() : Field<T>() {};
    // Methods
void allocate_memory(bool not_empty, int arr_sz) override {
  if (not_empty) {
    class_eval[0] = new double[arr_sz];
    class_eval[1] = new double[arr_sz];
    class_eval[2] = new double[arr_sz];
    }
  else {
    class_eval = {0, 0, 0};
    };
  }

void free_memory(bool not_empty) override {
  if (not_empty) {
    delete class_eval[0];
    delete class_eval[1];
    delete class_eval[2];
  };
}


public:

  ~RegularVectorField() {free_memory(has_grid);};

  // Constructors
  RegularVectorField () : Field<std::array<double, 3>, std::array<double*, 3>>() {
      size_t ar_sz = array_size();
      allocate_memory(has_grid, ar_sz);
      };

  RegularVectorField (std::array<int, 3>  shape, std::array<double, 3>  zeropoint, std::array<double, 3>  grid_increment) : Field<std::array<double, 3>, std::array<double*, 3>>(shape, zeropoint, grid_increment) {
      size_t ar_sz = array_size();
      allocate_memory(has_grid, ar_sz);
      };

  RegularVectorField (std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : Field<std::array<double, 3>, std::array<double*, 3>>(grid_x, grid_y, grid_z) {
      size_t ar_sz = array_size();
      allocate_memory(has_grid, ar_sz);
      };


  // Fields

  const int ndim = 3;
  // Methods

  std::array<double*, 3> on_grid(int seed = 0) {
    if (not has_grid) {
      throw GridException();
    }
    if (has_grid) { 
      if (regular_grid) {
        evaluate_function_on_grid(class_eval, shape, zeropoint, increment, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
        return class_eval;
      }
      else {
        evaluate_function_on_grid(class_eval, grid_x, grid_y, grid_z, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
        return class_eval;
      }
    }
    else
      return class_eval;

  }

  std::array<double*, 3> on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, int seed = 0) {
    std::array<double*, 3> function_eval;
    for (int i = 0; i<3; ++i) {
      function_eval[i] = new double[grid_x.size()*grid_y.size()*grid_z.size()];
    }
    std::cout << "The ref to b_grid when created 1 " << &function_eval << " \n\n";
    evaluate_function_on_grid(function_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    std::cout << "The ref to b_grid when created 2" << &function_eval << " \n\n";
    return function_eval;
  }

 std::array<double*, 3> on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, int seed = 0) {
    std::array<double*, 3> function_eval;
    for (int i = 0; i<3; ++i) {
      function_eval[i] = new double[grid_x.size()*grid_y.size()*grid_z.size()];
    }
    evaluate_function_on_grid(function_eval, grid_shape, grid_zeropoint, grid_increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return function_eval;
  }

};

#endif
