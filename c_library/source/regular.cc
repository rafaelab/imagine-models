#include <cmath>
#include <iostream>
#include "ImagineModels/MagneticField.h"
#include "ImagineModels/ThermalElectronField.h"

 std::vector<std::vector<std::vector<std::vector<double>>>> RegularMagneticField::_evaluate_grid(
    std::vector<double> grid_x, const std::vector<double> grid_y, const std::vector<double> grid_z,
    std::function<std::vector<double>(std::vector<double> )> ev_at_pos ) const {
    int size_x = grid_x.size();
    int size_y = grid_y.size();
    int size_z = grid_z.size();
//    if (size != grid_y.size()) {
//        throw std::length_error(std::string("size of grid y (" + std::to_string(grid_y.size())) +
//        std::string(") is not equal to size of grid x (" + std::to_string(size) + std::string(")")));}
//    if (size != grid_z.size()) {
//        throw std::length_error(std::string("size of grid z (" + std::to_string(grid_z.size())) +
//        std::string(") is not equal to size of grid x (" + std::to_string(size) + std::string(")")));}
    std::vector<std::vector<std::vector<std::vector<double>>>> b(size_x , std::vector<std::vector<std::vector<double>>> (size_y, std::vector<std::vector<double>> (size_z, std::vector<double> (3))));
    for (int i=0; i < size_x; i++) {
        for (int j=0; j < size_y; j++) {
            for (int k=0; k < size_z; k++) {
                std::vector<double> grid_at_pos = {{grid_x[i], grid_y[j], grid_z[k]}};
                b[i][j][k] = ev_at_pos(grid_at_pos);
            }
        }
    }
    return b;
  }


std::vector<std::vector<std::vector<double>>> RegularThermalElectronField::_evaluate_grid(
     std::vector<double> grid_x, const std::vector<double> grid_y, const std::vector<double> grid_z,
     std::function<double(std::vector<double> )> ev_at_pos ) const {
     int size_x = grid_x.size();
     int size_y = grid_y.size();
     int size_z = grid_z.size();
//     if (size != grid_y.size()) {
//         throw std::length_error(std::string("size of grid y (" + std::to_string(grid_y.size())) +
//         std::string(") is not equal to size of grid x (" + std::to_string(size) + std::string(")")));}
//     if (size != grid_z.size()) {
//         throw std::length_error(std::string("size of grid z (" + std::to_string(grid_z.size())) +
//         std::string(") is not equal to size of grid x (" + std::to_string(size) + std::string(")")));}
     std::vector<std::vector<std::vector<double>>> x(size_x , std::vector<std::vector<double>> (size_y, std::vector<double> (size_z)));
     for (int i=0; i < size_x; i++) {
         for (int j=0; j < size_y; j++) {
             for (int k=0; k < size_z; k++) {
                 std::vector<double> grid_at_pos = {{grid_x[i], grid_y[j], grid_z[k]}};
                 x[i][j][k] = ev_at_pos(grid_at_pos);
             }
         }
     }
     return x;
   }
