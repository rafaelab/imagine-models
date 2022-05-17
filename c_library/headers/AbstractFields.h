#include <vector>
#include <cassert>
#include <stdexcept>
#include <iostream>

#include "exceptions.h"


template<typename G>
class Field {
protected:
    bool regular;
    bool vector_valued;
public:
  // constructors
  Field() = default;
  virtual ~Field() {};

  Field(const G &grid_x, const G &grid_y, const G &grid_z):
      gx(grid_x), gy(grid_y), gz(grid_z), size(grid_x.size()) {
    assert(size == grid_y.size());
    assert(size == grid_z.size());
    };

  // Fields
  G &gx; // x coordinates
  G &gy; // y coordinates
  G &gz; // z coordinates
  int size;



};

template<typename G>
class SumField : protected Field<G>{
protected:
public:
  SumField(const Field<G> &summand1, const Field<G> &summand2) :
  Field<G>(summand1.gx, summand1.gy, summand1.gz),
    summand1(summand1), summand2(summand2)  {
      // Consistency checks
      assert(summand1.size == summand2.size);
      assert(summand1.gx == summand2.gx);
      assert(summand1.gy == summand2.gy);
      assert(summand1.gz == summand2.gz);
      if (summand1.regular && summand2.regular) {
        regular = true;
      } else {
        regular = false;
      }
      if (summand1.vector_valued && summand2.vector_valued) {
        vector_valued = true
      } else if {
        if (!summand1.vector_valued && !summand2.vector_valued) {
          vector_valued = false
      }
      else {
        throw FieldException();
      }
    }

  Field<G> &summand1;
  Field<G> &summand2;


  double *evaluate_model(const double x, const double y, const double z) {
    if constexpr (regular) {
      return *(*summand1.evaluate_model(x, y, z) + *summand2.evaluate_model(x, y, z));
    } else {
      throw NotImplementedException();
    }   }



};

template<typename G>
class RegularField : protected Field<G> {
protected:
  bool regular = true;
public:
  // constructors
  using RegularField<G> :: Field<G>;
  // methods

  RegularField operator+(const RegularField& f) {
         SumField sum(*this, f);
         return sum;
       };

  virtual double *evaluate_model(const double &x, const double &y, const double &z) = 0;

  virtual std::vector<double> getField(const std::vector<double> &pos_vec) {
    // This is the interface function to CRPRopa
    double *evm = evaluate_model(pos_vec[0], pos_vec[1], pos_vec[2]);
    return std::vector<double>(evm, evm + 3);
  };

  std::vector<double> grid_to_scalar_field(const G &ggx, const G &ggy, const G &ggz, const int size,
  std::function<double*(double, double, double)> eval) const{
     std::vector<double> b(size*size*size);
     for (int i=0; i < size; i++) {
         int m = i*size;
         for (int j=0; j < size; j++) {
            int n = j*size;
            for (int k=0; k < size; k++) {
              b[m + n + k] = *eval(ggx.at(i), ggy.at(j), ggz.at(k));
         }   }   }
     return b;
     }

   std::vector<double> grid_to_vector_field(const G &ggx, const G &ggy, const G &ggz, const int size
   std::function<double*(double, double, double)> eval) {
      std::vector<double> b(size*size*size*3);
      for (int i=0; i < size; i++) {
          int m = i*size;
          for (int j=0; j < size; j++) {
              int n = i*size;
              for (int k=0; k < size; k++) {
                  int o = k*size;
                  double v = *eval(ggx.at(i), ggy.at(j), ggz.at(k));
                  for (int l=0; l < 3; l++) {
                      b[m + n + o + l] = v[l];
          }   }   }
      return b;
      }


   std::vector<double> evaluate_model_on_grid(const G &grid_x, const G &grid_y, const G &grid_z) {
     int siz = grid_x.size();
     assert(siz == grid_y.size());
     assert(siz == grid_z.size());
     if constexpr (vector_valued) {
       std::vector<double> b = grid_to_vector_field(grid_x, grid_y, grid_z, siz,
       evaluate_model);
     } else {
       std::vector<double> b = grid_to_scalar_field(grid_x, grid_y, grid_z, siz,
       evaluate_model);
     }
       return b;
   }

   std::vector<double> evaluate_model_on_grid() {
     if constexpr (vector_valued) {
       std::vector<double> b = grid_to_vector_field(gx, gy, gz, size,
       evaluate_model);
     } else {
       std::vector<double> b = grid_to_scalar_field(gx, gy, gz, size,
       evaluate_model);
     }
     return b;
   }
};
