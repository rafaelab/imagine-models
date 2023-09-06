#include <functional>
#include <cmath>

//#include "param.h"
#include "Field.h"
#include "RegularField.h"

#include "param.h"


struct HelixParams : Params {
    //parameter_names = {"ampx", "ampy"};

    number ampx = 0; // all differentiable parameters are real
    number ampy = 0;
    number ampz = 0;

    double rmax = 20.;
    double rmin = 1.;
    };



class HelixMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
        bool has_derivative = true;
        vector  _at_position(const double &x, const double &y, const double &z, const HelixParams &p) const;

    public:
        using RegularVectorField :: RegularVectorField;

        HelixParams param;

        vector at_position(const double &x, const double &y, const double &z) const {
            return _at_position(x, y, z, this->param);
        }
        

        #if autodiff_FOUND
            Eigen::MatrixXd _derivative(const double &x, const double &y, const double &z,  HelixParams &p) {
                vector out;
                Eigen::MatrixXd deriv = ad::jacobian([&](auto _x, auto _y, auto _z, auto _p) {return this->_at_position(_x, _y, _z, _p);}, ad::wrt(p.ampx, p.ampy, p.ampz), ad::at(x, y, z, p), out);  
                return deriv;
            }

            Eigen::MatrixXd derivative(const double &x, const double &y, const double &z) {
                return _derivative(x, y, z, this->param);
            }
        #endif        
 };
