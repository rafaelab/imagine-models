#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "Field.h"
#include "RandomField.h"

class ESRandomField : public RandomVectorField {
  protected:
    bool DEBUG = false;
  public:
    using RandomVectorField :: RandomVectorField;

    double r0 = 8.5;
    double z0 = 1.5;
    std::array<double, 3> observer{8.5, 0, 0};


    void _on_grid(std::array<double*, 3> grid_eval, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) override;


    double spatial_profile(const double &x, const double &y, const double &z) const; 

};
