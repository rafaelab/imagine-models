#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>
#include <array>


inline vector Cyl2Cart(double phi, vector invec){
	vector outvec{{0., 0., 0.}};
	double cosphi = std::cos(phi);
	double sinphi = std::sin(phi);
	auto inv = cosphi * invec[0];
	outvec[0] = cosphi * invec[0] - sinphi * invec[1];
	outvec[1] = sinphi * invec[0] + cosphi * invec[1];
	outvec[2] = invec[2];
	return outvec;
}// Cyl2Cart


inline void addVector(vector vec, vector vec2add){
	vec[0] += vec2add[0];
	vec[1] += vec2add[1];
	vec[2] += vec2add[2];
}// addVector

#endif