#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "RegularModels.h"

void test_parameter_update() {
    UniformMagneticField umf = UniformMagneticField();
    assert (umf.bx == 0.);
    assert (umf.by == 0.);
    assert (umf.bz == 0.);

    std::array<double, 3> zeros{0., 0., 0.};
    
    assert (umf.at_position(2.4, 2.1, -.2) == zeros);
    
    umf.bx = -3.2; 
    
    assert (umf.bx == -3.2); 
    std::array<double, 3> updated{-3.2, 0., 0.};
    
    assert (umf.at_position(2.4, 2.1, -.2) == updated);

}

int main() {
    test_parameter_update();
}