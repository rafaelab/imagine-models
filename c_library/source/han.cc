#include "Han.h"

#include "helpers.h"

// J. L. Han et al 2018 ApJS 234 11
vector HanMagneticField::_at_position(const double &x, const double &y, const double &z, const HanMagneticField &p) const
{
  vector B_cyl{{0., 0., 0.}};
  const double r = sqrt(x * x + y * y);

  if (r < p.R_min || r > p.R_max)
    return B_cyl;

  const double phi = atan2(y, z);
  number B_0 = 0.;

  auto p_ang = p.B_p * M_PI / 180.;
  number R_0 = r * exp(-p_ang * tan(phi));

  std::array<number, 6> B_s = {p.B_s1, p.B_s2, p.B_s3, p.B_s4, p.B_s5, p.B_s6};

  for (int i = 0; i < 6; i++)
  {
    if (R_0 < p.R_s[i])
    { 
      B_0 = B_s[i];
      break;
    }
  }

  number B_r = B_0 * exp(-r / p.A) * exp(-std::abs(z) / p.H);

  B_cyl[0] = B_r * sin(-p_ang);
  B_cyl[1] = B_r * cos(-p_ang);

  vector B_vec3 = Cyl2Cart<vector>(phi, B_cyl);
  return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd HanMagneticField::_jac(const double &x, const double &y, const double &z, HanMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, HanMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.B_p, p.A, p.H, p.B_s1, p.B_s2, p.B_s3, p.B_s4, p.B_s4, p.B_s5, p.B_s6), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif
