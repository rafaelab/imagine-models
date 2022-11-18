#include <vector>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <algorithm>

#include "exceptions.h"


template<typename G, typename T>
class RegularField {
protected:
    // Fields
    // constructors
    RegularField() {};

    RegularField(const RegularField<G, T>& f) {};
    // methods

public:

  // Fields

  // methods
  virtual T evaluate_model(const double &x, const double &y, const double &z) const = 0;

  std::vector<double> getField(const std::vector<double> &pos_vec) {
    // This is the interface function to CRPRopa
    T evm = evaluate_model(pos_vec[0], pos_vec[1], pos_vec[2]);
    return evm;
  }

    // Evaluate scalar valued functions
  std::vector<double> evaluate_function_on_grid(const G &ggx, const G &ggy, const G &ggz, const int size,
  std::function<double(double, double, double)> eval) const{
     std::vector<double> b(size*size*size);
     for (int i=0; i < size; i++) {
         int m = i*size;
         for (int j=0; j < size; j++) {
             int n = j*size;
             for (int k=0; k < size; k++) {
                 b[m + n + k] = eval(ggx.at(i), ggy.at(j), ggz.at(k));
         }   }   }
     return b;
     }

     // Evaluate vector valued functions
   std::vector<double> evaluate_function_on_grid(const G &ggx, const G &ggy, const G &ggz,
     const int size, std::function<std::vector<double>(double, double, double)> eval) {
      std::vector<double> b(size*size*size*3);
      std::vector<double>::iterator bend = b.begin();
      for (int i=0; i < size; i++) {
          for (int j=0; j < size; j++) {
              for (int k=0; k < size; k++) {
                  std::vector<double> v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
                  bend = std::copy(v.begin(), v.end() + 3, bend);
          }   }   }
      return b;
    }


   std::vector<double> evaluate_model_on_grid(const G &grid_x, const G &grid_y, const G &grid_z) {
     int siz = grid_x.size();
     assert(siz == grid_y.size());
     assert(siz == grid_z.size());
     std::vector<double> b = evaluate_function_on_grid(grid_x, grid_y, grid_z, siz,
       evaluate_model);
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

  T evaluate_model(const double& x, const double& y, const double& z) const {

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
      std::vector<double> sum1 = summand1.evaluate_model(x, y, z);
      std::vector<double> sum2 = summand2.evaluate_model(x, y, z);
      std::vector<double> sum(3);
      for (size_t i; i<3; i++) {
        sum[i] = sum1[i] + sum2[i];
      }
      return sum;
      }

  double sum_scalar_models(const double& x, const double& y, const double& z) const {
      return summand1.evaluate_model(x, y, z) + summand2.evaluate_model(x, y, z);
    }

};
