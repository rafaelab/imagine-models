#include <functional>
#include <cmath>

//#include "param.h"
#include "Field.h"
#include "RegularField.h"
#include "autodiff.hh"
#include "param.h"

#if defined autodiff_FOUND
    #include <autodiff/forward/real.hpp>
    #include <autodiff/forward/dual.hpp>
    #include <autodiff/forward/real/eigen.hpp>
    namespace ad = autodiff;
    typedef ad::real number;
    typedef ad::VectorXreal vector;
#else
    typedef double number;  // only used for differentiable numbers! 
    typedef std::array<double, 3> vector;
#endif


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
    public:
        using RegularVectorField :: RegularVectorField;

        HelixParams param;
        

        std::array<double, 3>  at_position(const double &x, const double &y, const double &z) const;

        #if defined autodiff_FOUND

            vector  _at_position(const double &x, const double &y, const double &z, const HelixParams &p) const;

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
