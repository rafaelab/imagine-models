#include "../../ImagineModels/src/jf12.cc"
#include <cassert>
#include <iostream>

int main() {
  JF12MagneticField jf12;
  std::vector<double> test_pos{{1., 2., 3.}};
  std::vector<double> b_val(3, 1.e-12);
  std::vector<double> jf12_val = jf12.evaluate_at_pos(test_pos);
  double jf12_b_arm_1_val = jf12.b_arm_1;
  std::cout << jf12.b_arm_1 << std::endl;
  assert ((b_val ==  jf12_val));
  return 0;
}