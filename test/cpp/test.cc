#include "../../ImagineModels/c/headers/MagneticField.h"
#include <cassert>
#include <iostream>
#include <vector>

int main() {
  JF12MagneticField jf12;
  std::vector<double> test_pos{{6., 2., 0.}};
  std::vector<double> b_val = {1.e-12, 1.e-12, 1.e-12};
  std::vector<double> jf12_val = jf12.evaluate_at_pos(test_pos);
  double jf12_b_arm_1_val = jf12.b_arm_1;
  std::cout << "b_arm_1: " << jf12.b_arm_1 << std::endl;
  std::cout << "jf12_val at ";
  for (int l = 0; l < test_pos.size(); l++) {
                    std::cout << test_pos[l] << " kpc, ";
                    }
  std::cout << ":\n";
  for (int l = 0; l < jf12_val.size(); l++) {
                    std::cout << jf12_val[l]<< ", ";
                    }
  std::cout << "\n";
  jf12.b_arm_1 = 40.;
  std::cout << "updated b_arm_1: " << jf12.b_arm_1 << std::endl;
  std::cout << "updated jf12_val at ";
  for (int l = 0; l < test_pos.size(); l++) {
                    std::cout << test_pos[l] << " kpc, ";
                    }
  std::cout << ":\n";
  for (int l = 0; l < jf12_val.size(); l++) {
                    std::cout << jf12_val[l] << ", ";
                    }
  std::cout << "\n";
  const std::vector<double> grid_x {{2,4,6,8}};
  const std::vector<double> grid_y {{4,6,8,10}};
  const std::vector<double> grid_z {{-2,-4,-6,-8}};
  return 0;
}