#include <vector>
#include <stdexcept>

// magnetic field base class
class MagneticField {
public:
  MagneticField() = default;
//  MagneticField(const MagneticField &) = delete;
//  MagneticField(MagneticField &&) = delete;
//  MagneticField &operator=(const MagneticField &) = delete;
//  MagneticField &operator=(MagneticField &&) = delete;
  virtual ~MagneticField() = default;
};

// regular magnetic field
class RegularMagneticField : public MagneticField {
public:
  RegularMagneticField() = default;
//  RegularMagneticField(const RegularMagneticField &) = delete;
//  RegularMagneticField(RegularMagneticField &&) = delete;
//  RegularMagneticField &operator=(const RegularMagneticField &) = delete;
//  RegularMagneticField &operator=(RegularMagneticField &&) = delete;
  virtual ~RegularMagneticField() = default;

  virtual std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const {std::vector<double> _b(3); return _b;}

  std::vector<std::vector<std::vector<std::vector<double>>>> evaluate_at_grid(const std::vector<double> grid_x,
                                       const std::vector<double> grid_y,
                                       const std::vector<double> grid_z) const {
    int size = sizeof(grid_x)/sizeof(grid_x[0]);
    if (size != sizeof(grid_y)/sizeof(grid_y[0])) {
        throw std::length_error("size of grid y is not equal to size of grid x ");}
    if (size != sizeof(grid_z)/sizeof(grid_z[0])) {
        throw std::length_error("size of grid z is not equal to size of grid x ");}
    std::vector<std::vector<std::vector<std::vector<double>>>> b(size , std::vector<std::vector<std::vector<double>>> (size, std::vector<std::vector<double>> (size, std::vector<double> (3))));
    for (int i = 0; i < size; i++) {
        for (int j=0; j < size; j++) {
            for (int k=0; k < size; k++) {
                std::vector<double> grid_at_pos = {{grid_x[i], grid_y[j], grid_z[k]}};
                b[i][j][k] = evaluate_at_pos(grid_at_pos);
            }
        }
    }
    return b;
  }
 };
