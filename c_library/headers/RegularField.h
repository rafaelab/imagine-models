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
void allocate_memory(std::array<int, 3> shp, bool do_allocation, bool call_from_init) {
  if (do_allocation) {
    size_t arr_sz = array_size();
    grid_eval = new double[arr_sz];
    }
  else {
    grid_eval = 0;
    };
  }

void free_memory(bool do_deallocation) {
  if (do_deallocation) {
    delete grid_eval;
  };
}

public:
  // -----CONSTRUCTORS-----

  ~RegularScalarField() {free_memory(true);};
  
  RegularScalarField () : Field<double, double*>() {
      allocate_memory({0, 0, 0}, false, true);
      };

  RegularScalarField (std::array<int, 3>  shape, std::array<double, 3>  zeropoint, std::array<double, 3>  increment) : Field<double, double*>(shape, zeropoint, increment) {
      allocate_memory(shape, true, true);
      };

  RegularScalarField (std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : Field<double, double*>(grid_x, grid_y, grid_z) {
      allocate_memory(shape, true, true);
      };

  // Fields
  const int ndim = 1;
  // Methods

  double* on_grid(int seed = 0) {
    if (not initialized_with_grid) {
      throw GridException();
    }
    if (regular_grid) {
      evaluate_function_on_grid(grid_eval, shape, zeropoint, increment, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
      return grid_eval;
    }
    else {
      evaluate_function_on_grid(grid_eval, grid_x, grid_y, grid_z, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
      return grid_eval;
    }
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
void allocate_memory(std::array<int, 3> shp, bool do_allocation, bool call_from_init) override {
  if (do_allocation) {
    size_t arr_sz = array_size();
    grid_eval[0] = new double[arr_sz];
    grid_eval[1] = new double[arr_sz];
    grid_eval[2] = new double[arr_sz];
    }
  else {
    grid_eval = {0, 0, 0};
    };
  }

void free_memory(bool do_deallocation) override {
  if (do_deallocation) {
    delete grid_eval[0];
    delete grid_eval[1];
    delete grid_eval[2];
  };
}


public:

  ~RegularVectorField() {free_memory(true);};

  // Constructors
  RegularVectorField () : Field<std::array<double, 3>, std::array<double*, 3>>() {
      allocate_memory({0, 0, 0}, false, true);
      };

  RegularVectorField (std::array<int, 3>  shape, std::array<double, 3>  zeropoint, std::array<double, 3>  grid_increment) : Field<std::array<double, 3>, std::array<double*, 3>>(shape, zeropoint, grid_increment) {
      allocate_memory(shape, true, true);
      };

  RegularVectorField (std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : Field<std::array<double, 3>, std::array<double*, 3>>(grid_x, grid_y, grid_z) {
      allocate_memory(shape, true, true);
      };


  // Fields

  const int ndim = 3;
  // Methods

  std::array<double*, 3> on_grid(int seed = 0) {
    if (not initialized_with_grid) {
      throw GridException();
    }

    if (regular_grid) {
      evaluate_function_on_grid(grid_eval, shape, zeropoint, increment, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
      return grid_eval;
    }
    else {
      evaluate_function_on_grid(grid_eval, grid_x, grid_y, grid_z, [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
      return grid_eval;
    }

  }

  std::array<double*, 3> on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, int seed = 0) {
    std::array<double*, 3> function_eval;
    for (int i = 0; i<3; ++i) {
      function_eval[i] = new double[grid_x.size()*grid_y.size()*grid_z.size()];
    }
    evaluate_function_on_grid(function_eval, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
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
