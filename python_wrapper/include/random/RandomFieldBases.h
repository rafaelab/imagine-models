#ifndef RANDOMFIELDBASES_H
#define RANDOMFIELDBASES_H

#include <pybind11/pybind11.h>

#include "../random_trampoline.h"
#include "../array_converters.h"

namespace py = pybind11;
using namespace pybind11::literals;

void RandomFieldBases(py::module_ &m) {
    // Random Vector Base Class
    py::class_<RandomVectorField, RandomField<vector, std::array<double*, 3>>, PyRandomVectorField>(m, "RandomVectorField")
      .def(py::init<>())
      .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

      .def("on_grid", [](RandomVectorField &self, std::array<int, 3> &grid_shape,  std::array<double, 3>  &grid_reference_point, std::array<double, 3>  &grid_increment, int seed)  {
          std::array<double*, 3> f = self.on_grid(grid_shape, grid_reference_point, grid_increment, seed);
          size_t sx = grid_shape[0]; 
          size_t sy = grid_shape[1];
          size_t sz = grid_shape[2];
          
          auto lis = from_pointer_array_to_list_pyarray(f, sx, sy, sz, true);
          return lis;},
          py::kw_only(), py::arg("shape").noconvert(), py::arg("reference_point"), py::arg("increment"),  "seed"_a, 
          py::return_value_policy::take_ownership)

        .def("on_grid", [](RandomVectorField &self, int seed)  {
          std::array<double*, 3> f = self.on_grid(seed);
          size_t sx = self.shape[0];
          size_t sy = self.shape[1];
          size_t sz = self.shape[2]; // catches fftw zeropad (uneven)

          auto arr = from_pointer_array_to_list_pyarray(f, sx, sy, sz, true);

          return arr;}, 
          "seed"_a, 
          py::return_value_policy::take_ownership);

// Random Scalar Base Class
    py::class_<RandomScalarField, RandomField<number, double*>, PyRandomScalarField>(m, "RandomScalarField")
      .def(py::init<>())
      .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

      .def("on_grid", [](RandomScalarField &self, std::array<int, 3> &grid_shape,  std::array<double, 3>  &grid_reference_point, std::array<double, 3>  &grid_increment, int seed)  {
          double* f = self.on_grid(grid_shape, grid_reference_point, grid_increment, seed);
          size_t sx = grid_shape[0];
          size_t sy = grid_shape[1];
          size_t sz = grid_shape[2];

          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz, true);

          return arr;},
          py::kw_only(), py::arg("shape").noconvert(), py::arg("reference_point"), py::arg("increment"),  "seed"_a, 
          py::return_value_policy::take_ownership)

     .def("on_grid", [](RandomScalarField &self, int seed)  {
          double* f = self.on_grid(seed);
          size_t sx = self.shape[0];
          size_t sy = self.shape[1];
          size_t sz = self.shape[2]; // catches fftw zeropad (uneven)
          
          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz, true);
          
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;}, 
          "seed"_a, 
          py::return_value_policy::take_ownership);
      
}

#endif