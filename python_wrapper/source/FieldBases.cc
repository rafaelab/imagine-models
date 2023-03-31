void FieldBases(py::module_ &) {
    
    py::class_<Field<std::array<double, 3>, std::array<double*, 3>>,  PyVectorFieldBase>(m, "VectorFieldBase");

    py::class_<Field<double, double*>,  PyScalarFieldBase>(m, "ScalarFieldBase");

    py::class_<RandomField<std::array<double, 3>, std::array<double*, 3>>,  PyVectorRandomFieldBase>(m, "VectorRandomFieldBase");

    py::class_<RandomField<double, double*>,  PyScalarRandomFieldBase>(m, "ScalarRandomFieldBase");

}