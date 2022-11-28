
#include "../../c_library/headers/hamunits.h"
#include "../../c_library/headers/AbstractFields.h"
#include "../../c_library/headers/RegularJF12.h"
#include "../../c_library/headers/Jaffe.h"
#include "../../c_library/headers/Helix.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

void print_pos(std::map <std::string, std::vector<double>> pd,
               std::map <std::string, std::shared_ptr<RegularField<std::vector<double>, std::vector<double>>>> md) {

      auto model_iter = md.begin();


      while (model_iter != md.end()) {
        auto position_iter = pd.begin();
        std::cout << "The nodel " << model_iter->first << " is evaluated \n\n";
        while (position_iter != pd.end()) {
          std::vector<double> mval = (*(model_iter->second)).getField(position_iter->second);

          std::cout << "Position: ";
          for (size_t l = 0; l < (position_iter->second).size(); l++) {
                            std::cout << (position_iter->second)[l] << " kpc  ";
                          }
          std::cout << "\n";
          std::cout << "Evaluation: ";
          for (size_t l = 0; l < mval.size(); l++) {
                            std::cout << mval[l] << "  ";
                            }
          std::cout << "\n\n";

          ++position_iter;
          }
        std::cout << "\n";
        ++model_iter;
      }
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

  std::map <std::string, std::vector<double>> position_dict;


  // Define some positions in Galactic cartesian coordinates (units are kpc)
  position_dict["origin"] =  {0., 0., 0.};
  position_dict["pos_on_x_axis"] = {1.7, 0., 0.};
  position_dict["pos_on_y_axis"] = {0., -2.1, 0.};
  position_dict["pos_on_z_axis"] = {0., 0., 1.35};
  position_dict["pos_on_xy_plane"] = {1.9, -2.9, 0.};
  position_dict["pos_on_yz_plane"] = {0., 2.1, -1.6};
  position_dict["pos_on_xz_plane"] = {-2.7, 0., 1.25};

  // Define a grid in Galactic cartesian coordinates (units are kpc)
  const std::vector<double> grid_x {{2., 4., 0., 1., .4}};
  const std::vector<double> grid_y {{4., 6., 0.1, 0., .2}};
  const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 1.}};


  // using pointer here since RegularField is abstract
  std::map <std::string, std::shared_ptr<RegularField<std::vector<double>, std::vector<double>>>> model_dict;

  model_dict["Jansson Farrar"] = std::shared_ptr<JF12MagneticField<std::vector<double>>> (new JF12MagneticField<std::vector<double>>());
  model_dict["Jaffe"] = std::shared_ptr<JaffeMagneticField<std::vector<double>>> (new JaffeMagneticField<std::vector<double>>());
  model_dict["Helix"] = std::shared_ptr<HelixMagneticField<std::vector<double>>> (new HelixMagneticField<std::vector<double>>());


  print_pos(position_dict, model_dict);

  //std::vector<double> jf12_grid = jf12_1.evaluate_model_on_grid(grid_x, grid_y, grid_z);

  //std::vector<std::vector<std::vector<double>>> ymw_grid = ymw16.evaluate_grid(grid_x, grid_y, grid_z);

  //print_ev_grid(jf12_grid, grid_x, grid_y, grid_z);


//  std::cout << "\n b_arm_1: " << jf12.b_arm_1 << std::endl;



  //jf12.b_arm_1 = 40.;

  //std::cout << "\n updated b_arm_1: " << jf12.b_arm_1 << std::endl;

  // Evaluate model again
  //std::vector<double> jf12_val2 = jf12.getField(position_dict.at("pos_on_xy_plane"));
  return 0;
}
