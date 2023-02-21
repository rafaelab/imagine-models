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


template<typename T>
class RegularField : public Field<T>  {
protected:
    // Fields

    // Constructors
    using Field<T> :: Field;
    RegularField() : Field<T>() {};
    // Methods

public:
  // Constructors
  virtual ~RegularField() = default;
  // Fields

  // Methods
  virtual T at_position(const double &x, const double &y, const double &z) const = 0;

  on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z) {
    std::vector<int> siz{(int)grid_x.size(), grid_y.size(), grid_z.size()};
    std::vector<double> b;
    this->evaluate_function_on_grid(b, siz, grid_x, grid_y, grid_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return b;
  }

  on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) {
    double b;
    this->evaluate_function_on_grid(b, n, zeropoint, increment,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
    return b;
  }

};


class RegularScalarField: public RegularField<double> {};
class RegularVectorField: public RegularField<std::vector<double>> {};

#endif
