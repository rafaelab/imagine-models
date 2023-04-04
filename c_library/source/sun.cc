#include <cmath>
#include "../headers/hamunits.h"
#include "../headers/Sun.h"

#include "../headers/helpers.h"

//Sun et al. A&A V.477 2008 ASS+RING model magnetic field
std::array<double, 3>  Sun2008MagneticField::at_position(const double &x, const double &y, const double &z) const {

double r = sqrt(x*x + y*y);
double phi = atan2(y, z);
  
  // first we set D2
  // ------------------------------------------------------------
  double D2;
  if(r > 7.5)
    {
      D2 = 1.;
    }
  else if(r <= 7.5 && r > 6.)
    {
      D2 = -1.;
    }
  else if(r <= 6. && r > 5.)
    {
      D2 = 1.;
    }
  else if(r <= 5.)
    {
      D2 =- 1.;
    }
  else {std::cerr << " Error! r in Sun field is:" << r << "kpc " << std::endl; exit(1);}
  // ------------------------------------------------------------

  // now we set D1
  // ------------------------------------------------------------
  double D1;
  if(r > b_Rc)
    {
      D1 = b_B0 * exp(-((r - b_Rsun) / b_R0) - (std::abs(z) / b_z0));
    }
  else if(r <= b_Rc)
    {
      D1 = b_Bc;
    }
  else {std::cerr << " Error! r in Sun field is:" << r << "kpc " << std::endl; exit(1);}
  // ------------------------------------------------------------
  
  double p_ang = b_pitch_deg * M_PI / 180.;
  std::array<double, 3>  B_cyl{D1 * D2 * std::sin(p_ang),
		                           -D1 * D2 * std::cos(p_ang),
		                           0.};

  // [ORIGINAL HAMMURABI COMMENT]  Taking into account the halo field
  double halo_field;

  // [ORIGINAL HAMMURABI COMMENT]  for better overview
  double bH_z1;
  if (std::abs(z) < bH_z0) {
    bH_z1 = bH_z1a; 
    }
  else { 
    bH_z1 = bH_z1b;
  }
  double hf_piece1 = (bH_z1 * bH_z1) / (bH_z1 * bH_z1 + (std::abs(z) - bH_z0) * (std::abs(z) - bH_z0));
  double hf_piece2 = exp(-(r-bH_R0) / (bH_R0));

  halo_field = bH_B0 * hf_piece1 * (r / bH_R0) * hf_piece2;

  // [ORIGINAL HAMMURABI COMMENT] Flip north.  Not sure how Sun did this.  This is his code with no
  // flip though the paper says it's flipped, but without this mod,
  // there is no antisymmetry across the disk.  However, it doesn't seem to work.  
  if (z > 0 && r >= 5.) { 
    halo_field*=-1.;
    }

  B_cyl[1] += halo_field;

  std::array<double, 3> B_vec3;

  B_vec3 = Cyl2Cart(phi, B_cyl);

  return B_vec3; 
};