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
    // constructors
    Field() {};

    Field(const Field<G, T>& f): regular(f.regular) {
      };
    // methods
    void set_regular(bool reg) {
      regular = reg;
    }

public:

  // Fields
  bool regular;
  // methods
  virtual T evaluate_model(const double &x, const double &y, const double &z) const = 0;

  std::vector<double> getField(const std::vector<double> &pos_vec) {
    // This is the interface function to CRPRopa
    T evm = evaluate_model(pos_vec[0], pos_vec[1], pos_vec[2]);
    return evm;
  }

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


template<typename G>
class VectorField : public Field<G, std::vector<double>> {
protected:
  VectorField() : Field<G, std::vector<double>>() {};
public:
  // constructors
  using Field<G, std::vector<double>> :: Field;
};


template<typename G>
class SumVectorField : public VectorField<G>{
public:
  SumVectorField(const VectorField<G>& summand1, const VectorField<G>& summand2) :
  VectorField<G>(summand1),
    summand1(summand1), summand2(summand2)  {
      // Consistency checks
      if (summand1.regular && summand2.regular) {
        VectorField<G>::set_regular(true);
      } else {
        VectorField<G>::set_regular(false);
      }
      }

  const VectorField<G> &summand1;
  const VectorField<G> &summand2;

  std::vector<double> evaluate_model(const double& x, const double& y,
    const double& z) const {
      std::vector<double> sum1 = summand1.evaluate_model(x, y, z);
      std::vector<double> sum2 = summand2.evaluate_model(x, y, z);
      std::vector<double> sum(3);
      for (size_t i; i<3; i++) {
        sum[i] = sum1[i] + sum2[i];
      }
      return sum;
      }

};



template<typename G>
class ScalarField : public Field<G, double> {
protected:
public:
    // constructors
    using Field<G, double> :: Field;
    ScalarField();
    virtual ~ScalarField() {};

};


template<typename G>
class SumScalarField : public ScalarField<G>{
public:
  SumScalarField(const ScalarField<G>& summand1, const ScalarField<G>& summand2) :
  ScalarField<G>(summand1),
    summand1(summand1), summand2(summand2)  {
      // Consistency checks
      if (summand1.regular && summand2.regular) {
        ScalarField<G>::set_regular(true);
      } else {
        ScalarField<G>::set_regular(false);
      }
      }

  const ScalarField<G> &summand1;
  const ScalarField<G> &summand2;

  double evaluate_model(const double& x, const double& y, const double& z) const {
      return summand1.evaluate_model(x, y, z) + summand2.evaluate_model(x, y, z);
    }
};
