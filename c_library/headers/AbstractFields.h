#include <vector>
#include <cassert>
#include <stdexcept>
#include <iostream>


template<typename G>
class Field {
protected:
    bool regular;
    bool vector_valued;
public:
  // constructors
  Field() = default;
  virtual ~Field() {}

  Field(const G &grid_x, const G &grid_y, const G &grid_y):
      gx(grid_x), gy(grid_y),size(grid_x.size()) {
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
class SumField : Field<G>{
private:
public:
  SumField(const Field<G> &summand1, const Field<G> &summand2) :
    gx(summand1.gx), gy(summand1.gy), gz(summand1.gz), size(summand1.size),
    summand1(summand1), summand2(summand2) {
      assert(size == summand2.size);
      assert(gx == summand2.gx);
      assert(gy == summand2.gy);
      assert(gz == summand2.gz);
    }

  Field &summand1;
  Field &summand2;


  double ev_at_pos(const double x, const double y, const double z) {
    return summand1.ev_at_pos(x, y, z) + summand2.ev_at_pos(x, y, z);
  }



};

template<typename G>
class RegularField : Field {
protected:
  bool regular = true;
public:
  // constructors
  using RegularField  :: Field;
  // methods

  double *evaluate_model(const double x, const double y, const double z) = 0;

  std::vector<double> grid_to_scalar_field(const G &ggx, const G &ggy, const G &ggz, const int size,
  std::function< *double(double, double, double)> eval) const{
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
   std::function< *double(double, double, double)> eval) {
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
          }   }  }
      return b;
      }


   std::vector<double> evaluate_model_at_positions(const G &grid_x, const G &grid_y, const G &grid_z) {
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

   std::vector<double> evaluate_positions() {
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
