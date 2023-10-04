#include <cmath>

#include "TinyakovTkachev.h"
#include "helpers.h"

vector TTMagneticField::_at_position(const double &x, const double &y, const double &z, const TTMagneticField &p) const
{

    double r = sqrt(x * x + y * y);

    vector B_vec3{{0, 0, 0}};
    if (r <= b_r_max)
    {
        return B_vec3;
    }

    double phi = atan2(y, x);

    auto beta = 1. / tan(p.b_p);
    if (abs(p.b_p) < 1.e-16)
    {
        beta = 1.;
    }

    auto phase = (beta * log(1. + p.b_d / p.b_Rsun)) - M_PI / 2.;
    //  double epsilon0=(b5_Rsun+b5_d)*exp(-(M_PI/2.)*tan(b5_p)); <-- hammurabi comment

    double sign;
    if (z < 0)
    {
        sign = 1.;
    }
    else
    {
        sign = -1.;
    }
    auto f_z = sign * exp(-(abs(z) / p.b_z0));

    // there is a factor 1/cos(phase) difference between
    // the original TT and Kachelriess. <-- hammurabi comment
    // double b_r=b_b0*(b_Rsun/(r));
    auto b_r = p.b_b0 * (p.b_Rsun / (r * cos(phase)));

    // if(r<b5_r_min) {b_r=b5_b0*(b5_Rsun/(b5_r_min));} <-- hammurabi comment
    if (r < b_r_min)
    {
        b_r = p.b_b0 * (p.b_Rsun / (b_r_min * cos(phase)));
    }

    auto B_r_phi = b_r * cos(phi - beta * log(r / p.b_Rsun) + phase);
    //  if(r<1.e-26){B_r_phi = b_r*std::cos(phi+phase);} <-- hammurabi comment

    // B-field in cylindrical coordinates: <-- hammurabi comment
    vector B_cyl{{B_r_phi * sin(p.b_p) * f_z,
                  B_r_phi * cos(p.b_p) * f_z,
                  0.}};

    B_vec3 = Cyl2Cart<vector>(phi, B_cyl);

    return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd TTMagneticField::_jac(const double &x, const double &y, const double &z, TTMagneticField &p) const
{
    vector out;
    Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, TTMagneticField &_p)
                                          { return _p._at_position(_x, _y, _z, _p); },
                                          ad::wrt(p.b_Rsun, p.b_b0, p.b_d, p.b_z0, p.b_p), ad::at(x, y, z, p), out);
    return _filter_diff(_deriv);
}

#endif