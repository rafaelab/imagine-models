#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>
#include <array>


template <typename ARRAY>
inline void Cyl2Cart(double phi, ARRAY invec, ARRAY outvec){
	double cosphi = cos(phi);
	double sinphi = sin(phi);
	outvec[0] = cosphi * invec[0] - sinphi * invec[1];
	outvec[1] = sinphi * invec[0] + cosphi * invec[1];
	outvec[2] = invec[2];
}// Cyl2Cart

#endif