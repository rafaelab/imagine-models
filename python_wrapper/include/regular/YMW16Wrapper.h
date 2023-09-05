#include <pybind11/pybind11.h>

#include "YMW.h"

namespace py = pybind11;
using namespace pybind11::literals;

void YMW(py::module_ &m) {

    py::class_<YMWParams>(m, "YMWParams")
        .def(py::init<>())
        .def_readwrite("r_warp", &YMWParams::r_warp)
        .def_readwrite("r0", &YMWParams::r0)
        .def_readwrite("t0_gamma_w", &YMWParams::t0_gamma_w)
        // Thick disc
        .def_readwrite("t1_ad", &YMWParams::t1_ad)
        .def_readwrite("t1_bd", &YMWParams::t1_bd)
        .def_readwrite("t1_n1", &YMWParams::t1_n1)
        .def_readwrite("t1_h1", &YMWParams::t1_h1)
        // Thin disc
        .def_readwrite("t2_a2", &YMWParams::t2_a2)
        .def_readwrite("t2_b2", &YMWParams::t2_b2)
        .def_readwrite("t2_n2", &YMWParams::t2_n2)
        .def_readwrite("t2_k2", &YMWParams::t2_k2)
        // Spiral arms
        .def_readwrite("t3_b2s", &YMWParams::t3_b2s)
        .def_readwrite("t3_ka", &YMWParams::t3_ka)
        .def_readwrite("t3_aa", &YMWParams::t3_aa)
        .def_readwrite("t3_ncn", &YMWParams::t3_ncn)
        .def_readwrite("t3_wcn", &YMWParams::t3_wcn)
        .def_readwrite("t3_thetacn", &YMWParams::t3_thetacn)
        .def_readwrite("t3_nsg", &YMWParams::t3_nsg)
        .def_readwrite("t3_wsg", &YMWParams::t3_wsg)
        .def_readwrite("t3_thetasg", &YMWParams::t3_thetasg)
        .def_readwrite("t3_rmin", &YMWParams::t3_rmin)
        .def_readwrite("t3_phimin", &YMWParams::t3_phimin)
        .def_readwrite("t3_tpitch", &YMWParams::t3_tpitch)
        .def_readwrite("t3_cpitch", &YMWParams::t3_cpitch)
        .def_readwrite("t3_narm", &YMWParams::t3_narm)
        .def_readwrite("t3_warm", &YMWParams::t3_warm)
        // Galactic center
        .def_readwrite("t4_ngc", &YMWParams::t4_ngc)
        .def_readwrite("t4_agc", &YMWParams::t4_agc)
        .def_readwrite("t4_hgc", &YMWParams::t4_hgc)
        // Gum
        .def_readwrite("t5_kgn", &YMWParams::t5_kgn)
        .def_readwrite("t5_ngn", &YMWParams::t5_ngn)
        .def_readwrite("t5_wgn", &YMWParams::t5_wgn)
        .def_readwrite("t5_agn", &YMWParams::t5_agn)
        // Local bubble
        .def_readwrite("t6_j_lb", &YMWParams::t6_j_lb)
        .def_readwrite("t6_nlb1", &YMWParams::t6_nlb1)
        .def_readwrite("t6_detlb1", &YMWParams::t6_detlb1)
        .def_readwrite("t6_wlb1", &YMWParams::t6_wlb1)
        .def_readwrite("t6_hlb1", &YMWParams::t6_hlb1)
        .def_readwrite("t6_thetalb1", &YMWParams::t6_thetalb1)
        .def_readwrite("t6_nlb2", &YMWParams::t6_nlb2)
        .def_readwrite("t6_detlb2", &YMWParams::t6_detlb2)
        .def_readwrite("t6_wlb2", &YMWParams::t6_wlb2)
        .def_readwrite("t6_hlb2", &YMWParams::t6_hlb2)
        .def_readwrite("t6_thetalb2", &YMWParams::t6_thetalb2)
        // Loop
        .def_readwrite("t7_nli", &YMWParams::t7_nli)
        .def_readwrite("t7_rli", &YMWParams::t7_rli)
        .def_readwrite("t7_wli", &YMWParams::t7_wli)
        .def_readwrite("t7_detthetali", &YMWParams::t7_detthetali)
        .def_readwrite("t7_thetali", &YMWParams::t7_thetali);

    py::class_<YMW16, RegularScalarField>(m, "YMW16")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())

        .def("at_position", &YMW16::at_position, "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
        
        .def_readwrite("param", &YMW16::param);

}