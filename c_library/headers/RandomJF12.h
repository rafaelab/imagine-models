#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "Field.h"
#include "RandomField.h"

class JF12RandomField : public RandomVectorField {
  protected:
    bool DEBUG = false;
  public:
    using RandomVectorField :: RandomVectorField;

    JF12RandomField() : RandomVectorField() {};
    ~JF12RandomField() {};

    double b0_1 = 10.81; // uG
    double b0_2 = 6.96; // uG
    double b0_3 = 9.59; // uG
    double b0_4 = 6.96; // uG
    double b0_5 = 1.96; // uG
    double b0_6 = 16.34; // uG
    double b0_7 = 37.29; // uG
    double b0_8 = 10.35; // uG
    double b0_int = 7.63; // uG
    double z0_spiral = 0.61; // kpc

    double b0_halo = 4.68; // uG
    double r0_halo = 10.97; // kpc
    double z0_halo = 2.84; // kpc

    double Rmax = 20.;
    double rho_GC = 1.;

    double spectral_amplitude = 1.; 
    double spectral_offset = 1.; 
    double spectral_slope = 2.;

    void _on_grid(std::array<double*, 3> grid_eval, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed);

    double calculate_fourier_sigma(const double &abs_k) const override;

    double spatial_profile(const double &x, const double &y, const double &z) const; 
};
