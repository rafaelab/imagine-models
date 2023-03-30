#include <cmath>
#include <array>



void Cyl2Cart( double phi, std::array<double, 3> &invec, std::array<double, 3> &outvec){
	outvec[0]=0.;outvec[1]=0.;outvec[2]=0.;
	std::array<std::array<double, 3>, 3> cyl_unit_vec_array;
	cyl_unit_vec_array[0][0]=cos(phi);
	cyl_unit_vec_array[1][0]=sin(phi);
	cyl_unit_vec_array[2][0]=0.;
	cyl_unit_vec_array[0][1]=-sin(phi);
	cyl_unit_vec_array[1][1]=cos(phi);
	cyl_unit_vec_array[2][1]=0.;
	cyl_unit_vec_array[0][2]=0.;
	cyl_unit_vec_array[1][2]=0.;
	cyl_unit_vec_array[2][2]=1.;

	for (int n=0;n<3;n++){
		for (int m=0;m<3;m++) {
			outvec[n]=outvec[n]+cyl_unit_vec_array[n][m]*invec[m];
		}
	}
}// Cyl2Cart

void Cyl2Cart2( double phi, double invec[3], double outvec[3]){
	double cosphi=cos(phi);
	double sinphi=sin(phi);
	outvec[0]=cosphi*invec[0]-sinphi*invec[1];
	outvec[1]=sinphi*invec[0]+cosphi*invec[1];
	outvec[2]=invec[2];
}// Cyl2Cart
