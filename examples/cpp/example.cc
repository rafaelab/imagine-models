#include "../../c_library/headers/MagneticField.h"
#include <cassert>
#include <iostream>
#include <vector>

void print_ev_pos(JF12MagneticField m_to_eval, std::vector<double> tp) {
      std::vector<double> jf12_val = m_to_eval.evaluate_at_pos(tp);
      std::cout << "jf12_val at ";
      for (int l = 0; l < tp.size(); l++) {
                        std::cout << tp[l] << " kpc, ";
                        }
      std::cout << ":\n";
      for (int l = 0; l < jf12_val.size(); l++) {
                        std::cout << jf12_val[l]<< ", ";
                        }
      std::cout << "\n";
      }

int main() {
  JF12MagneticField jf12;
  std::vector<double> test_pos{{6., 2., 0.}};
  std::vector<double> b_val = {1.e-12, 1.e-12, 1.e-12};
  const std::vector<double> grid_x {{2.,4.,6.,8., 0., 0., 0., 0.}};
  const std::vector<double> grid_y {{4.,6.,8.,2., 0., 0., 0., 0.}};
  const std::vector<double> grid_z {{-.2,-.4, .6, 0, -.1, .8, .2, 0}};
  std::vector<std::vector<std::vector<std::vector<double>>>>  b_grid = jf12.evaluate_grid(grid_x, grid_y, grid_z);
  std::cout<< "b_grid: "  <<std::endl;
  for (int i = 0; i < b_grid.size(); i++) {
        std::cout<< "b_grid x size: " << b_grid.size() <<std::endl;
        for (int j = 0; j < b_grid[i].size(); j++) {
            std::cout<< "b_grid y size: " << b_grid[i].size() <<std::endl;
            for (int k = 0; k < b_grid[i][j].size(); k++) {
                std::cout<< "b_grid z size: " << b_grid[i][j].size() <<std::endl;
                std::cout << "x: "<< grid_x[i] << ", " << "y: "<< grid_y[j]<< ", "<< "z: "<< grid_z[k] << ", ";
                for (int l = 0; l < b_grid[i][j][k].size(); l++) {
                    std::cout << b_grid[i][j][k][l] << " ";
                    }
                std::cout<<"\n";
                }
            }
        }
  std::cout << "\n b_arm_1: " << jf12.b_arm_1 << std::endl;
  print_ev_pos(jf12, test_pos);

  // Change some parameters
  jf12.b_arm_1 = 40.;
  jf12.b_arm_2 = 40.;
  jf12.b_arm_3 = 40.;
  jf12.b_arm_4 = 40.;
  jf12.b_arm_5 = 40.;
  jf12.b_arm_6 = 40.;
  jf12.b_arm_7 = 40.;
  std::cout << "b_arm_1: " << jf12.b_arm_1 << std::endl;
  print_ev_pos(jf12, test_pos);
  return 0;
}
