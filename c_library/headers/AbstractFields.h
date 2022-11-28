#include <vector>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fftw3.h>

// #include "exceptions.h"



template<typename G, typename T>
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
  virtual std::vector<double> on_grid(const G &grid_x, const G &grid_y, const G &grid_z) const = 0;

  // This is the interface function to CRPRopa
  std::vector<double> getField(const std::vector<double> &pos_vec) const {
    T evm = at_position(pos_vec[0], pos_vec[1], pos_vec[2]);
    return evm;
  }
    // Evaluate scalar valued functions
  std::vector<double> evaluate_function_on_grid(const G &ggx, const G &ggy, const G &ggz,
    const int size_x, const int size_y, const int size_z,
  std::function<double(double, double, double)> eval) const {
     std::vector<double> b(size_x*size_y*size_z);
     for (int i=0; i < size_x; i++) {
         int m = i*size_x;
         for (int j=0; j < size_y; j++) {
             int n = j*size_y;
             for (int k=0; k < size_z; k++) {
                 b[m + n + k] = eval(ggx.at(i), ggy.at(j), ggz.at(k));
         }   }   }
     return b;
     }

     // Evaluate vector valued functions
   std::vector<double> evaluate_function_on_grid(const G &ggx, const G &ggy, const G &ggz,
     const int size_x, const int size_y, const int size_z,
     std::function<std::vector<double>(double, double, double)> eval) const {
      std::vector<double> b;
      std::vector<double>::iterator bend = b.begin();
      for (int i=0; i < size_x; i++) {
          for (int j=0; j < size_y; j++) {
              for (int k=0; k < size_z; k++) {
                  std::vector<double> v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
                  b.insert(b.end(), v.begin(), v.end());
          }   }   }
      return b;
    }
};


template<typename G, typename T>
class RegularField : public Field<G, T>  {
protected:
    // Fields

    // Constructors
    using Field<G, T> :: Field;
    RegularField() : Field<G, T>() {};
    // Methods

public:
  // Constructors
  virtual ~RegularField() = default;
  // Fields

  // Methods
  virtual T at_position(const double &x, const double &y, const double &z) const = 0;

  std::vector<double> on_grid(const G &grid_x, const G &grid_y, const G &grid_z) const {
     int siz_x = grid_x.size();
     int siz_y = grid_y.size();
     int siz_z = grid_z.size();

     std::vector<double> b = RegularField::evaluate_function_on_grid(grid_x, grid_y, grid_z, siz_x, siz_y, siz_z,
      [this](double xx, double yy, double zz) {return at_position(xx, yy, zz);});
  return b;
 }

};


template<typename G, typename T>
class SumRegularField : public RegularField<G, T>{
public:
  SumRegularField(const RegularField<G, T>& summand1, const RegularField<G, T>& summand2) :
  RegularField<G, T>(summand1),
    summand1(summand1), summand2(summand2)  {
      }

  const RegularField<G, T> &summand1;
  const RegularField<G, T> &summand2;

  T at_position(const double& x, const double& y, const double& z) const {

  T sum;
  if constexpr (std::is_same<T, double>::value) {
    sum = sum_scalar_models(x, y, z);
    }
  else
    {
    sum =  sum_vector_models(x, y, z);
  }
  return sum;
}

  std::vector<double> sum_vector_models(const double& x, const double& y,
    const double& z) const {
      std::vector<double> sum1 = summand1.at_position(x, y, z);
      std::vector<double> sum2 = summand2.at_position(x, y, z);
      std::vector<double> sum(3);

      std::transform(sum1.begin(), sum1.end(), sum2.begin(), sum.begin(), std::plus<double>());
      return sum;
      }

  double sum_scalar_models(const double& x, const double& y, const double& z) const {
      return summand1.at_position(x, y, z) + summand2.at_position(x, y, z);
    }

};

template<typename G, typename T>
class RandomField  : public Field<G, T>  {
protected:
    // Fields
    // constructors
    using Field<G, T> :: Field;

    RandomField() : Field<G, T>() {};
    ~RandomField() {};

    // methods

public:

  // Fields

  // methods
  virtual T at_position(const double &x, const double &y, const double &z) const = 0;

  std::vector<double> on_grid(const G &grid_x, const G &grid_y, const G &grid_z) = 0;


  virtual spatial_profile(const double &x, const double &y, const double &z) const = 0;

  // this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
  // original author: https://github.com/gioacchinowang
  double spectrum(const double &k, const double rms, const double k0, const double k1, const double a0, const double a1) const {
  const double p0 = rms*rms;
  const double unit = 1. / (4 * cgs::pi * k * k);   // units fixing, wave vector in 1/kpc units
  // power laws
  const double band1{double(k < k1)};
  const double band2{double(k > k1) * double(k < k0)};
  const double band3{double(k > k0)};
  const double P =
      band1 *
          std::pow(k0 / k1, a1) *
          std::pow(k / k1, 6.0) +
      band2 / std::pow(k / k0, a1) +
      band3 / std::pow(k / k0, a0);
  return P * p0 * unit;
}

// this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
// original author: https://github.com/gioacchinowang
std::vector<double> divergence_cleaner(const Hamvec<3, ham_float> &k,
                       const Hamvec<3, ham_float> &b) const {

  double k_length = 0;
  double b_length = 0;
  double b_dot_k = 0;
  for (int i = 0; i < 3; ++i) {
    k_length += static_cast<double>(k[i] * k[i]);
    b_length += static_cast<double>(b[i] * b[i]);
    b_dot_k += static_cast<double>(b[i] * k[i]);
    }
  if (k_length == 0 or b_length == 0) {
    return b;
    }
  k_length = std::sqrt(k_length);

  const double inv_k_mod = 1. / k_length;
  // multiply \sqrt(3) for preserving spectral power statistically
  return std::vector<double>{
      1.73205081 * (b[0] - k[0] * b_dot_k * inv_k_mod),
      1.73205081 * (b[1] - k[1] * b_dot_k * inv_k_mod),
      1.73205081 * (b[2] - k[2] * b_dot_k * inv_k_mod)};
}


};
