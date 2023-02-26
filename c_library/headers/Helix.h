#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

class HelixMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;

        double ampx = 0.;
        double ampy = 0.;
        double ampz = 0.;
        double rmax = 3.;
        double rmin = 0.;

        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const {
          const double r{sqrt(x*x + y*y)}; // radius in cylindrical coordinates
          const double phi{atan2(y, x) + M_PI / 2.0}; // azimuthal angle in cylindrical coordinates#
          std::array<double, 3> b{0.0, 0.0, 0.0};
          if ((r > rmin) && (r < rmax)) {
            b = std::array<double, 3> {ampx * std::cos(phi), ampy * std::sin(phi), ampz};
            }
          return b;
        };
        //std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override;
        /*
        std::vector<double> _dampx_at_pos(const std::vector<double> &pos) const;
        std::vector<double> _dampy_at_pos(const std::vector<double> &pos) const;
        std::vector<double> _dampz_at_pos(const std::vector<double> &pos) const;

        std::vector<std::vector<std::vector<std::vector<double>>>> dampx_grid(const std::vector<double> grid_x,
                                          const std::vector<double> grid_y,
                                          const std::vector<double> grid_z) const {
                                          return _evaluate_grid(grid_x, grid_y, grid_z, [this](std::vector<double> p)
                                          {return _dampx_at_pos(p);});
                                          };
        std::vector<std::vector<std::vector<std::vector<double>>>> dampy_grid(const std::vector<double> grid_x,
                                          const std::vector<double> grid_y,
                                          const std::vector<double> grid_z) const {
                                          return _evaluate_grid(grid_x, grid_y, grid_z, [this](std::vector<double> p)
                                          {return _dampy_at_pos(p);});
                                        };
        std::vector<std::vector<std::vector<std::vector<double>>>> dampz_grid(const std::vector<double> grid_x,
                                          const std::vector<double> grid_y,
                                          const std::vector<double> grid_z) const {
                                          return _evaluate_grid(grid_x, grid_y, grid_z, [this](std::vector<double> p)
                                          {return _dampz_at_pos(p);});
                                                                          };
*/
 };
