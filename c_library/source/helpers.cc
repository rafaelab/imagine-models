#include <cmath>
#include <array>





void Cyl2Cart( double phi, double invec[3], double outvec[3]){
	double cosphi=cos(phi);
	double sinphi=sin(phi);
	outvec[0]=cosphi*invec[0]-sinphi*invec[1];
	outvec[1]=sinphi*invec[0]+cosphi*invec[1];
	outvec[2]=invec[2];
}// Cyl2Cart
