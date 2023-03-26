#include <ImagineModels/RegularModels.h>
#include <ImagineModels/RandomModels.h>

#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>



  std::map <std::string, std::array<double, 3>> position_dict;


  // Define some positions in Galactic cartesian coordinates (units are kpc)
  position_dict["origin"] =  {0., 0., 0.};
  position_dict["pos_on_x_axis"] = {1.7, 0., 0.};
  position_dict["pos_on_y_axis"] = {0., -2.1, 0.};
  position_dict["pos_on_z_axis"] = {0., 0., 1.35};
  position_dict["pos_on_xy_plane"] = {1.9, -2.9, 0.};
  position_dict["pos_on_yz_plane"] = {0., 2.1, -1.6};
  position_dict["pos_on_xz_plane"] = {-2.7, 0., 1.25};

  // Define a irregular grid in Galactic cartesian coordinates (units are kpc)
  const std::vector<double> grid_x {{2., 4., 0., 1., .4}};
  const std::vector<double> grid_y {{4., 6., 0.1, 0., .2}};
  const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 1.}};


  // Define a regular grid in Galactic cartesian coordinates (units are kpc)
  const std::array<int, 3> shape {{4, 3, 2}};
  const std::array<double, 3> increment {{2.1, 0.3, 1.}};
  const std::array<double, 3> zeropoint {{-4., 0.1, -0.3}};

  // Dictionaries
  // using pointer here since RegularField is abstract
  std::map <std::string, std::shared_ptr<RegularVectorField>> reg_mods;

  reg_mods["Jansson Farrar regular"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField());
  reg_mods["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField());
  reg_mods["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField());


regular_models = ['JaffeMagneticField', 'HelixMagneticField', 'JF12RegularField',]

known_positions = {'JF12RegularField': [([0, 0, 0], [0, 0, 0]), # Galactic center
                                        ([.1, .3, .4], [0, 0, 0]), # Within inner boundary
                                          ([3, 20., .3], [0, 0, 0, ])], # outside outer boundary 
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

        