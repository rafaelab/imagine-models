#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "RegularModels.h"

void test_parameter_update() {
    UniformMagneticField umf = UniformMagneticField();
    assert (umf.param.bx == 0.);
    assert (umf.param.by == 0.);
    assert (umf.param.bz == 0.);

    vector zeros{{0., 0., 0.}};
    
    assert (umf.at_position(2.4, 2.1, -.2) == zeros);
    
    umf.param.bx = -3.2; 
    
    assert (umf.param.bx == -3.2); 
    vector updated{{-3.2, 0., 0.}};
    
    assert (umf.at_position(2.4, 2.1, -.2) == updated);

}

int main() {
    test_parameter_update();
}