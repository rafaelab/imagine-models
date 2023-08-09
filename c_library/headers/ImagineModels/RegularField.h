#ifndef REGULARFIELD_H
#define REGULARFIELD_H


#include <vector>
#include <array>
#include <stdexcept>
#include <functional>
#include <iostream>
#include <algorithm>

#include "ImagineModels/exceptions.h"
#include "ImagineModels/Field.h"

class RegularScalarField : public Field<double, double*>  {
protected:
    // Fields

    // RegularField() : Field<T>() {};
    // Methods
double* allocate_memory(std::array<int, 3> shp) {
  size_t arr_sz = grid_size(shp);
  double* grid_eval = new double[arr_sz];
  return grid_eval;
  }

void free_memory(double* grid_eval) {
  delete grid_eval;
}

public:
  // -----CONSTRUCTORS-----

  ~RegularScalarField() {};
  
  RegularScalarField () : Field<double, double*>() {
      };

  RegularScalarField (std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  increment) : Field<double, double*>(shape, reference_point, increment) {
      };

  RegularScalarField (std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : Field<double, double*>(grid_x, grid_y, grid_z) {
      };

  // Fields
  const int ndim = 1;
  // Methods

  double* on_grid(int seed = 0) {
    if (not initialized_with_grid) {
      throw GridException();
    }
    double* grid_eval = allocate_memory(shape);
    if (regular_grid) {
      evaluate_function_on_grid(grid_eval, shape, reference_point, increment, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    }
    else {
      evaluate_function_on_grid(grid_eval, grid_x, grid_y, grid_z, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    }
    return grid_eval;
  }

  double* on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, const int seed = 0) {
    double *grid_eval =  allocate_memory(shape);
    //auto grid_eval = std::make_shared<double>();
    evaluate_function_on_grid(grid_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return grid_eval;
  }

  double* on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_reference_point, const std::array<double, 3> &grid_increment, const int seed = 0) {
    double *grid_eval =  allocate_memory(shape);
    evaluate_function_on_grid(grid_eval, grid_shape, grid_reference_point, grid_increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return grid_eval;
  }

};


class RegularVectorField : public Field<std::array<double, 3>, std::array<double*, 3>>  {
protected:
    // Fields

    // RegularField() : Field<T>() {};
    // Methods
std::array<double*, 3> allocate_memory(std::array<int, 3> shp) override {
  std::array<double*, 3> grid_eval;
  size_t arr_sz = grid_size(shp);
  grid_eval[0] = new double[arr_sz];
  grid_eval[1] = new double[arr_sz];
  grid_eval[2] = new double[arr_sz];
  return grid_eval;
  }

void free_memory(std::array<double*, 3> grid_eval) override {
  delete grid_eval[0];
  delete grid_eval[1];
  delete grid_eval[2];
}


public:

  ~RegularVectorField() {
    };

  // Constructors
  RegularVectorField () : Field<std::array<double, 3>, std::array<double*, 3>>() {
      };

  RegularVectorField (std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  grid_increment) : Field<std::array<double, 3>, std::array<double*, 3>>(shape, reference_point, grid_increment) {
      };

  RegularVectorField (std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : Field<std::array<double, 3>, std::array<double*, 3>>(grid_x, grid_y, grid_z) {
      };


  // Fields

  const int ndim = 3;
  // Methods

  std::array<double*, 3> on_grid(int seed = 0) {
    if (not initialized_with_grid) {
      throw GridException();
    }
    std::array<double*, 3> grid_eval = allocate_memory(shape);;
    if (regular_grid) {
      evaluate_function_on_grid(grid_eval, shape, reference_point, increment, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    }
    else {
      evaluate_function_on_grid(grid_eval, grid_x, grid_y, grid_z, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    }
  return grid_eval;
  }

  std::array<double*, 3> on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, int seed = 0) {
    std::array<int, 3> grid_shape = {(int)grid_x.size(), (int)grid_y.size(), (int)grid_z.size()};
    std::array<double*, 3> grid_eval = allocate_memory(grid_shape);
    evaluate_function_on_grid(grid_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return grid_eval;
  }

 std::array<double*, 3> on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_reference_point, const std::array<double, 3> &grid_increment, int seed = 0) {
    std::array<double*, 3> grid_eval = allocate_memory(grid_shape);
    evaluate_function_on_grid(grid_eval, grid_shape, grid_reference_point, grid_increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return grid_eval;
  }

};

#endif
