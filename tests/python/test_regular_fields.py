import ImagineModels as img

import pytest
import numpy as np


regular_models = ['JaffeMagneticField', 'HelixMagneticField', 'JF12RegularField', 'AxiSymmetricSpiral']

known_positions = {'JF12RegularField': [([0, 0, 0], [0, 0, 0]), # Galactic center
                                        ([.1, .3, .4], [0, 0, 0]), # Within inner boundary
                                          ([], [])], 
                   'JaffeMagneticField': [([0, 0, 0], [0, 0, 0]), # origin
                                          ([], [])], 
                   'AxiSymmetricSpiral': [([0, 0, 0], [0, 0, 0]), # origin
                                          ([], [])], 
                   'HelixMagneticField': [([0, 0, 0], [0, 0, 0]), # origin
                                          ([], [])]}


irregular_grid = {'grid_x': [-2.2, -1, 0, .1], 'grid_y': [-1, 0, 3.], 'grid_z': [-5, 0, .1]}

regular_grid = {'grid_shape': [2, 3, 4], 'increment': [.1, 3, .1], 'zeropoint': [-2., 3, .1]}

def test_at_position():
    for model_string in regular_models:
        mo = getattr(img, model_string) 
        mo_1 = mo()
        mo_2 = mo(**irregular_grid)
        mo_3 = mo(**regular_grid)
        
        for position, value in known_positions[model_string]:
            assert mo_1.at_position(*position) == value
            assert mo_2.at_position(*position) == value
            assert mo_3.at_position(*position) == value
            
            
def test_grid():
    # this function tests the consistency of the different 
    for model_string in regular_models:
        mo = getattr(img, model_string) 
        mo_1 = mo()
        mo_2 = mo(**irregular_grid)
        mo_3 = mo(**regular_grid)
        
        with pytest.raises(RuntimeError):
            mo_1.on_grid()
            
        assert mo_1.on_grid(**regular_grid) == mo_2.on_grid()
        assert mo_2.on_grid(**irregular_grid) == mo_3.on_grid()
            
        
def test_interface():
    raise NotImplementedError
            

def test_parameter_update():
    umf = img.UniformMagneticField()
    assert umf.bx == 0.
    assert umf.by == 0.
    assert umf.bz == 0.
    
    assert umf.at_position(2.4, 2.1, -.2) == [0., 0., 0]
    
    umf.bx = -3.2 
    
    assert umf.bx == -3.2 
    
    assert umf.at_position(2.4, 2.1, -.2) == [-3.2, 0., 0]

        