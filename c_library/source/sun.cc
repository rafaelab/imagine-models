#include <cmath>
#include "hamunits.h"
#include "Sun.h"

#include "helpers.h"

//Sun et al. A&A V.477 2008 ASS+RING model magnetic field
vector SunMagneticField::_at_position(const double &x, const double &y, const double &z, const SunParams &p) const {

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
  else //if(r <= 5.)
    {
      D2 =- 1.;
    }
  // ------------------------------------------------------------

  // now we set D1
  // ------------------------------------------------------------
  number D1;
  if(r >  p.b_Rc)
    {
      D1 =  p.b_B0 * exp(-((r -  p.b_Rsun) /  p.b_R0) - (abs(z) /  p.b_z0));
    }
  else //if(r <= b_Rc)
    {
      D1 =  p.b_Bc;
    }
  // ------------------------------------------------------------
  
  number p_ang = p.b_pitch_deg * M_PI / 180.;
  vector B_cyl{{ D1 * D2 * sin(p_ang), 
                -D1 * D2 * cos(p_ang), 
                 0. }};

  // [ORIGINAL HAMMURABI COMMENT]  Taking into account the halo field
  number halo_field;

  // [ORIGINAL HAMMURABI COMMENT]  for better overview
  number b3H_z1_actual;
  if (abs(z) <  p.bH_z0) {
    b3H_z1_actual =  p.bH_z1a;
    }
  else{b3H_z1_actual = p.bH_z1b;}
  number hf_piece1 = (b3H_z1_actual * b3H_z1_actual) / (b3H_z1_actual * b3H_z1_actual + (abs(z) - p.bH_z0) * (abs(z) - p.bH_z0));
  number hf_piece2 = exp(-(r - p.bH_R0) / (p.bH_R0));

  halo_field = p.bH_B0 * hf_piece1 * (r / p.bH_R0) * hf_piece2;

  // [ORIGINAL HAMMURABI COMMENT] Flip north.  Not sure how Sun did this.  This is his code with no
  // flip though the paper says it's flipped, but without this mod,
  // there is no antisymmetry across the disk.  However, it doesn't seem to work.  
  if (z > 0 && r >= 5.) { 
    halo_field*=-1.;
    }

  B_cyl[1] += halo_field;

  vector B_vec3;

  B_vec3 = Cyl2Cart(phi, B_cyl);

  return B_vec3; 
};