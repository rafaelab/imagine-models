#include <ImagineModels/RegularModels.h>
#include <ImagineModels/RandomModels.h>

#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>


int main() {
    std::map<std::string, std::map<std::array<double, 3>, std::array<double, 3>>> val_pos_map;

    // Define some positions in Galactic cartesian coordinates (units are kpc)

    std::map<std::array<double, 3>, std::array<double, 3>> jf_12_map;
    jf_12_map[{0., 0., 0.}] = {0., 0., 0.}; // Galactic center
    jf_12_map[{.1, .3, .4}] = {0., 0., 0.}; // Within inner boundary
    jf_12_map[{0., 0., 0.}] = {0., 0., 0.}; // outside outer boundary 


    std::map<std::array<double, 3>, std::array<double, 3>> jaffe_map;
    //jaffe_map[{0., 0., 0.}] = {0., 0., 0.}; // Galactic center

    std::map<std::array<double, 3>, std::array<double, 3>> helix_map;

    val_pos_map["JF12"] = jf_12_map;
    val_pos_map["Jaffe"] = jaffe_map;
    val_pos_map["Helix"] = helix_map;


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
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_empty_constructor;
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_regular_constructor;
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_irregular_constructor;

    models_w_empty_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField());
    models_w_empty_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField());
    models_w_empty_constructor["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField());

    models_w_regular_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(shape, increment, zeropoint));
    models_w_regular_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(shape, increment, zeropoint));
    models_w_regular_constructor["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField(shape, increment, zeropoint));

    models_w_empty_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(grid_x, grid_y, grid_z));
    models_w_empty_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(grid_x, grid_y, grid_z));
    models_w_empty_constructor["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField(grid_x, grid_y, grid_z));
}



void test_at_position(std::map<std::string, std::map<std::array<double, 3>, std::array<double, 3>>> val_pos_map,
                      std::map <std::string, std::shared_ptr<RegularVectorField>> model_dict
                      ) {
    auto model_iter = model_dict.begin();


    while (model_iter != model_dict.end()) {
        std::map<std::array<double, 3>, std::array<double, 3>> val_pos_map_this_model = val_pos_map[(model_iter->first)];
        auto val_pos_iter = val_pos_map_this_model.begin();
        while (val_pos_iter != val_pos_map_this_model.end()) {
          double x = (val_pos_iter->first)[0];
          double y = (val_pos_iter->first)[1];
          double z = (val_pos_iter->first)[2];

          std::array<double, 3> mval = (*(model_iter->second)).at_position(x, y, z);
          assert (mval != (val_pos_iter->second));
          ++val_pos_iter;
          }
        ++model_iter;
    }
}

void test_grid(std::map <std::string, std::shared_ptr<RegularVectorField>> models_no_grid, 
               std::map <std::string, std::shared_ptr<RegularVectorField>> models_regular_grid, 
               std::map <std::string, std::shared_ptr<RegularVectorField>> models_irregular_grid, 
               std::array<int, 3> shape, std::array<double, 3> increment, std::array<double, 3> zeropoint,
               std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z
              ) {
    auto model_iter = models_no_grid.begin();

    


    while (model_iter != models_no_grid.end()) { 
        std::array<double*, 3> eval_irregular = (*models_irregular_grid[model_iter->first]).on_grid(); 
        std::array<double*, 3> eval_regular = (*models_regular_grid[model_iter->first]).on_grid(); 
        std::array<double*, 3> eval_no_grid_regular = (*models_no_grid[model_iter->first]).on_grid(shape, increment, zeropoint); 
        std::array<double*, 3> eval_no_grid_irregular = (*models_no_grid[model_iter->first]).on_grid(grid_x, grid_y, grid_z); 
        size_t irreg_n = grid_x.size()*grid_y.size()*grid_z.size();
        size_t reg_n = shape[0]*shape[1]*shape[2];
        // std::array<double*, 3> eval_no_grid = (*models_no_grid[model_iter->first]).on_grid(); should raise an exception

        _check_array_equality_from_pointer(eval_irregular, eval_no_grid_irregular, irreg_n);
        _check_array_equality_from_pointer(eval_regular, eval_no_grid_regular, reg_n);
        ++model_iter;
    }

}

void _check_array_equality_from_pointer(std::array<double*, 3> a, std::array<double*, 3> b , size_t &n) {
    for (int d = 0; d==3; ++d) {
        std::vector<double>  arr_a(a[d], a[d] + n);
        std::vector<double>  arr_b(b[d], b[d] + n);
        assert (arr_a == arr_b); 
    }

}

            

            
        
//void test_interface() {
//}

            

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