#include <cmath>
#include <cassert>
#include <iostream>
#include "hamunits.h"
#include "RegularJF12.h"

vector JF12MagneticField::_at_position(const double &x, const double &y, const double &z, const JF12RegularParams &p) const {
  // define fixed parameters
      const double Rmax = 20;   // outer boundary of GMF
      const double rho_GC = 1.; // interior boundary of GMF

      // fixed disk parameters
      const double inc = 11.5; // inclination, in degrees
      const double rmin = 5.; // outer boundary of the molecular ring region
      const double rcent = 3. ; // inner boundary of the molecular ring region (field is
                         // zero within this region)
      const number f[8] = {
          0.130, 0.165, 0.094, 0.122,
          0.13,  0.118, 0.084, 0.156}; // fractions of circumference spanned by each
                                       // spiral, sums to unity
      const double rc_B[8] = {
          5.1,  6.3,  7.1, 8.3, 9.8,
          11.4, 12.7, 15.5}; // the radii where the spiral arm boundaries cross the
                             // negative x-axis

      const double r{sqrt(x * x + y * y)};
      const double rho{
          sqrt(x * x + y * y + z * z)};
      const double phi{atan2(y, x)};

      // define boundaries for where magnetic field is zero (outside of galaxy)
      if (r > Rmax || rho < rho_GC) {
        return vector{{0., 0., 0.}};
      }

      //------------------------------------------------------------------------------
      // DISK COMPONENT (8 spiral regions, 7 free parameters with 8th set to
      // conserve flux)
      // B0 set to 1 at r=5kpc
      const double B0 = (rmin / r); //
      // the logistic equation, to be multiplied to the toroidal halo field and
      // (1-zprofile) multiplied to the disk:
      const number zprofile{1. / (1 + exp(-2. / p.w_disk * (abs(z) - p.h_disk)))};

      // printf("%g, %g \n", z, zprofile);
      number B_cyl_disk[3] = {0, 0, 0}; // the disk field in cylindrical coordinates

      if ((r > rcent)) // disk field zero elsewhere
      {
        if (r < rmin) { // circular field in molecular ring
          B_cyl_disk[1] = B0 * p.b_ring * (1 - zprofile);
        } else {
          // use flux conservation to calculate the field strength in the 8th spiral
          // arm
          number bv_B[8] = {p.b_arm_1, p.b_arm_2, p.b_arm_3, p.b_arm_4,
                            p.b_arm_5, p.b_arm_6, p.b_arm_7, 0.};
          number b8 = 0.;

          for (int i = 0; i < 7; i++) {
            b8 -= f[i] * bv_B[i] / f[7];
          }
          bv_B[7] = b8;

          // iteratively figure out which spiral arm the current coordinates (r.phi)
          // correspond to
          number b_disk = 0.;
          double r_negx =
              r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi - M_PI));

          if (r_negx > rc_B[7]) {
            r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + M_PI));
          }
          if (r_negx > rc_B[7]) {
            r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + 3 * M_PI));
          }
          for (int i = 7; i >= 0; i--) {
            if (r_negx < rc_B[i]) {
              b_disk = bv_B[i];
            }
          } // "region 8,7,6,..,2"

          B_cyl_disk[0] = b_disk * B0 * sin(M_PI / 180. * inc) * (1 - zprofile);
          B_cyl_disk[1] = b_disk * B0 * cos(M_PI / 180. * inc) * (1 - zprofile);
        }
      }

      //-------------------------------------------------------------------------
      ////TOROIDAL HALO COMPONENT

      number b1, rh;
      number B_h = 0.;

      if (z >= 0) { // North
        b1 = p.Bn;
        rh = p.rn; // transition radius between inner-outer region
      } else {   // South
        b1 = p.Bs;
        rh = p.rs;
      }

      B_h = b1 * (1. - 1. / (1. + exp(-2. / p.wh * (r - rh)))) *
            exp(-(abs(z)) / (p.z0)); // vertical exponential fall-off
      const number B_cyl_h[3] = {0., B_h * zprofile, 0.};

      //------------------------------------------------------------------------
      // X- FIELD

      number Xtheta = 0.;
      number rp_X =  0.; // the mid-plane radius for the field line that pass through r
      number B_X = 0.;
      double r_sign = 1.; // +1 for north, -1 for south
      if (z < 0) {
        r_sign = -1.;
      }

      // dividing line between region with constant elevation angle, and the
      // interior:
      number rc_X = p.rpc_X + abs(z) / tan(p.Xtheta_const);

      if (r < rc_X) { // interior region, with varying elevation angle
        rp_X = r * p.rpc_X / rc_X;
        B_X = p.B0_X * pow(p.rpc_X / rc_X, 2.) * exp(-rp_X / p.r0_X);
        Xtheta = atan(abs(z) /
                      (r - rp_X)); // modified elevation angle in interior region
        // printf("Xtheta %g at z %g , r %g , rc_X %g \n",Xtheta, z,r,rc_X);
        if (z == 0.) {
          Xtheta = M_PI / 2.;
        }      // to avoid some NaN
      } else { // exterior region with constant elevation angle
        Xtheta = p.Xtheta_const;
        rp_X = r - abs(z) / tan(Xtheta);
        B_X = p.B0_X * rp_X / r * exp(-rp_X / p.r0_X);
      }

      // X-field in cylindrical coordinates
      number B_cyl_X[3] = {B_X * cos(Xtheta) * r_sign, 0., B_X * sin(Xtheta)};

      if (DEBUG) {

        std::cout << "DEBUG JF12: x: " <<  x <<  std::endl;
        std::cout << "DEBUG JF12: y: " <<  y <<  std::endl;
        std::cout << "DEBUG JF12: z: " <<  z <<  "\n" << std::endl;

        std::cout << "DEBUG JF12: r: " <<  r <<  std::endl;
        std::cout << "DEBUG JF12: rho: " <<  rho << std::endl;
        std::cout << "DEBUG JF12: phi: " <<  phi << "\n" << std::endl;

        std::cout << "DEBUG JF12: B_cyl_disk[0]: " <<  B_cyl_disk[0] << std::endl;
        std::cout << "DEBUG JF12: B_cyl_disk[1]: " <<  B_cyl_disk[1] << std::endl;
        std::cout << "DEBUG JF12: B_cyl_disk[2]: " <<  B_cyl_disk[2] <<  "\n" << std::endl;

        std::cout << "DEBUG JF12: B_cyl_h[0]: " <<  B_cyl_h[0] << std::endl;
        std::cout << "DEBUG JF12: B_cyl_h[1]: " <<  B_cyl_h[1] << std::endl;
        std::cout << "DEBUG JF12: B_cyl_h[2]: " <<  B_cyl_h[2] <<  "\n" << std::endl;

        std::cout << "DEBUG JF12: B_cyl_X[0]: " <<  B_cyl_X[0] << std::endl;
        std::cout << "DEBUG JF12: B_cyl_X[1]: " <<  B_cyl_X[1] << std::endl;
        std::cout << "DEBUG JF12: B_cyl_X[2]: " <<  B_cyl_X[2] <<  "\n" << std::endl;

      }


      // add fields together
      vector B_cyl{{0.0, 0.0, 0.0}};
      B_cyl[0] = B_cyl_disk[0] + B_cyl_h[0] + B_cyl_X[0];
      B_cyl[1] = B_cyl_disk[1] + B_cyl_h[1] + B_cyl_X[1];
      B_cyl[2] = B_cyl_disk[2] + B_cyl_h[2] + B_cyl_X[2];

      // convert field to cartesian coordinates
      vector B_cart{{0.0, 0.0, 0.0}};
      B_cart[0] = B_cyl[0] * cos(phi) - B_cyl[1] * sin(phi);
      B_cart[1] = B_cyl[0] * sin(phi) + B_cyl[1] * cos(phi);
      B_cart[2] = B_cyl[2];

      return B_cart;
  }
