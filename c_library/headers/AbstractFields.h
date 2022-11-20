#include <vector>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <algorithm>

#include "exceptions.h"



template<typename G, typename T>
class Field {
protected:
    // Fields
    // Constructors
    Field() {};

    Field(const Field<G, T>& f) {};
    // Methods

public:

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
      std::vector<double> b(size_x*size_y*size_z*3);
      std::vector<double>::iterator bend = b.begin();
      for (int i=0; i < size_x; i++) {
          for (int j=0; j < size_y; j++) {
              for (int k=0; k < size_z; k++) {
                  std::vector<double> v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
                  bend = std::copy(v.begin(), v.end() + 3, bend);
          }   }   }
      return b;
    }
};


template<typename G, typename T>
class RegularField : public Field<G, T>  {
protected:
    // Fields
    // constructors
    using Field<G, T> :: Field;

    RegularField() : Field<G, T>() {};
    virtual ~RegularField() {};

    // methods

public:

  // Fields

  // methods
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
      std::cout << "SumRegularField: sum1 " << sum1[0] << std::endl;
      std::cout << "SumRegularField: sum2 " << sum2[0] << std::endl;
      std::vector<double> sum(3);

      std::transform(sum1.begin(), sum1.end(), sum2.begin(), sum.begin(), std::plus<double>());
      //for (size_t i; i<3; i++) {
      //  sum[i] = sum1[i] + sum2[i];
      //}
      std::cout << "SumRegularField: sum " << sum[0] << std::endl;
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
    virtual ~RandomField() {};

    // methods

public:

  // Fields

  // methods
  virtual T at_position(const double &x, const double &y, const double &z) const = 0;

  std::vector<double> on_grid(const G &grid_x, const G &grid_y, const G &grid_z) = 0;

};
