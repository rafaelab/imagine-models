#include <cmath>
#include "hamunits.h"
#include "Helix.h"

vector HelixMagneticField::_at_position(const double &xx, const double &yy, const double &zz, const HelixMagneticField &p
) const {
  
  double phi = std::atan2(yy, xx); // azimuthal angle in cylindrical coordinates
  double r = std::sqrt(xx*xx + yy*yy); // radius in cylindrical coordinates
  vector b {{0.0, 0.0, 0.0}};
  if ((r > p.rmin) && (r < p.rmax)) {
    b[0] = std::cos(phi) * p.ampx; 
    b[1] = std::sin(phi) * p.ampy; 
    b[2] = p.ampz;
    }
  return b;
}

#if autodiff_FOUND

  Eigen::MatrixXd HelixMagneticField::_jac(const double &xx, const double &yy, const double &zz, HelixMagneticField &p) const {
    vector out;
    
    Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, HelixMagneticField &_p) {return _p._at_position(_x, _y, _z, _p);}, ad::wrt(p.ampx, p.ampy, p.ampz), ad::at(xx, yy, zz, p), out);  

    if (p.active_diff.size() != p.all_diff.size()) {
        std::vector<int> i_to_keep;
        
        for (std::string s : p.active_diff ) {
            if (auto search = p.all_diff.find(s); search != p.all_diff.end()) {
                int index = std::distance(p.all_diff.begin(), search);
                i_to_keep.push_back(index);
                }
        }
        return _deriv(Eigen::all, i_to_keep);
    }
    return _deriv;

};

#endif
