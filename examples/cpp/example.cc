#include "../../c_library/headers/MagneticField.h"
// #include "../../c_library/headers/ThermalElectronField.h"
#include <cassert>
#include <iostream>
#include <vector>

void print_ev_pos(std::vector<double> mval, std::vector<double> tp, const char* title) {

      std::cout << "Magnetic field at " << title << " (";
      for (size_t l = 0; l < tp.size(); l++) {
                        std::cout << tp[l] << " kpc, ";
                        }
      std::cout << "):\n";
      for (size_t l = 0; l < mval.size(); l++) {
                        std::cout << mval[l]<< ", ";
                        }
      std::cout << "\n";
      }

void print_ev_grid(std::vector<std::vector<std::vector<std::vector<double>>>>  b_grid,
  std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) {

  std::cout<< "b_grid: "  <<std::endl;
  for (size_t i = 0; i < b_grid.size(); i++) {
        std::cout<< "b_grid x size: " << b_grid.size() <<std::endl;
        for (size_t j = 0; j < b_grid[i].size(); j++) {
            std::cout<< "b_grid y size: " << b_grid[i].size() <<std::endl;
            for (size_t k = 0; k < b_grid[i][j].size(); k++) {
                std::cout<< "b_grid z size: " << b_grid[i][j].size() <<std::endl;
                std::cout << "x: "<< grid_x[i] << ", " << "y: "<< grid_y[j]<< ", "<< "z: "<< grid_z[k] << ", ";
                for (size_t l = 0; l < b_grid[i][j][k].size(); l++) {
                    std::cout << b_grid[i][j][k][l] << " ";
                    }
                std::cout<<"\n";
                }
            }
        }
    }

int main() {
  // initialize model
  JF12MagneticField jf12;

  // Define some positions in Galactic cartesian coordinates (units are kpc)
  std::vector<double> test_pos{{1., 2., 0.}};
  std::vector<double> origin{{0., 0., 0.}};
  std::vector<double> pos_on_x_axis{{.123, 0., 0.}};
  std::vector<double> pos_on_y_axis{{0., .432, 0.}};
  std::vector<double> pos_on_z_axis{{0., 0., 1.2}};
  std::vector<double> pos_on_xy_axis{{.123, -.213, 0.}};
  std::vector<double> pos_on_yz_axis{{0., -.232, -.3}};
  std::vector<double> pos_on_xz_axis{{-.24, 0., 1.2}};

  // evaluate the positions
  std::vector<double> jf12_test_pos = jf12.evaluate_at_pos(test_pos);
  print_ev_pos(jf12_test_pos, test_pos, "test_pos");
  std::vector<double> jf12_origin = jf12.evaluate_at_pos(origin);
  print_ev_pos(jf12_origin, origin, "origin");
  std::vector<double> jf12_pos_on_x_axis = jf12.evaluate_at_pos(pos_on_x_axis);
  print_ev_pos(jf12_pos_on_x_axis, pos_on_x_axis, "pos_on_x_axis");
  std::vector<double> jf12_pos_on_y_axis = jf12.evaluate_at_pos(pos_on_y_axis);
  print_ev_pos(jf12_pos_on_y_axis, pos_on_y_axis, "pos_on_y_axis");
  std::vector<double> jf12_pos_on_z_axis = jf12.evaluate_at_pos(pos_on_z_axis);
  print_ev_pos(jf12_pos_on_z_axis, pos_on_z_axis, "pos_on_z_axis");
  std::vector<double> jf12_pos_on_xy_axis = jf12.evaluate_at_pos(pos_on_xy_axis);
  print_ev_pos(jf12_pos_on_xy_axis, pos_on_xy_axis, "pos_on_xy_plane");
  std::vector<double> jf12_pos_on_yz_axis = jf12.evaluate_at_pos(pos_on_yz_axis);
  print_ev_pos(jf12_pos_on_yz_axis, pos_on_yz_axis, "pos_on_yz_plane");
  std::vector<double> jf12_pos_on_xz_axis = jf12.evaluate_at_pos(pos_on_xz_axis);
  print_ev_pos(jf12_pos_on_xz_axis, pos_on_xz_axis, "pos_on_xz_plane");

  // Define grid axes in Galactic cartesian coordinates (again, units are kpc)
  const std::vector<double> grid_x{2., 4., 0., 1.};
  const std::vector<double> grid_y{4., 6., 0.1, 0.};
  const std::vector<double> grid_z{-0.2, 0.8, 0.2, 0.};

  // evaluate the grid
  std::vector<std::vector<std::vector<std::vector<double>>>>  jf12_grid = jf12.evaluate_grid(grid_x, grid_y, grid_z);

  print_ev_grid(jf12_grid, grid_x, grid_y, grid_z);

  // print a parameter
  std::cout << "\n b_arm_1: " << jf12.b_arm_1 << std::endl;


  // Change the parameter
  jf12.b_arm_1 = 40.;

  std::cout << "\n updated b_arm_1: " << jf12.b_arm_1 << std::endl;

  // Evaluate moel again
  std::vector<double> jf12_val2 = jf12.evaluate_at_pos(test_pos);
  print_ev_pos(jf12_val2, test_pos, "updated jf12, test_pos");
  return 0;
}
