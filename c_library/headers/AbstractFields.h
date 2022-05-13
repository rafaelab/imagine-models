#include <vector>
#include <cassert>
#include <stdexcept>
#include <iostream>


template<typename G>
class Field {
private:
public:
  // Fields
  G gx; // x coordinates
  G gy; // y coordinates
  G gz; // z coordinates
  int size;

  // constructors
  Field() = default;
  virtual ~Field() {}

  Field(const G &grid_x, const G &grid_y, const G &grid_y):
      gx(grid_x), gy(grid_y),size(grid_x.size()) {
    assert(size == grid_y.size());
    assert(size == grid_z.size());
    };




  // methods

  std::vector<double> grid_to_scalar(const G &ggx, const G &ggy, const G &ggz, const int size) {
     std::vector<double> b(size*size*size);
     for (int i=0; i < size; i++) {
         int m = i*size;
         for (int j=0; j < size; j++) {
            int n = i*size;
            for (int l=0; l < size; l++) {
              b[m + n + l] = ev_at_pos(ggx.at(i), ggy.at(j), , ggz.at(j));
         }   }
     }
     return b;
     }


   std::vector<double> evaluate_grid(const G &grid_x, const G &grid_y, const G &grid_z) {
     int siz = grid_x.size();
     assert(siz == grid_y.size());
     assert(siz == grid_z.size());
     std::vector<double> b = grid_to_vector(grid_x, grid_y, grid_z, siz);
     return b;
   }

   std::vector<double> evaluate_grid() {
     std::vector<double> b = grid_to_vector(gx, gy, gz, size);
     return b;
   }



};
