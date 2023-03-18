#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "Field.h"
#include "RandomField.h"

class GaussianScalarField : public RandomScalarField {
  protected:
    bool DEBUG = false;
  public:
    using RandomScalarField :: RandomScalarField;

    void _on_grid(double* grid_eval, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) override;

    double spatial_profile(const double &x, const double &y, const double &z) const {
        return 1.;
    }; 

};
