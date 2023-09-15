#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>
#include <array>


inline vector Cyl2Cart(double phi, vector invec){
	vector outvec;
	double cosphi = std::cos(phi);
	double sinphi = std::sin(phi);
	outvec[0] = cosphi * invec[0] - sinphi * invec[1];
	outvec[1] = sinphi * invec[0] + cosphi * invec[1];
	outvec[2] = invec[2];
	return outvec;
}// Cyl2Cart

#endif