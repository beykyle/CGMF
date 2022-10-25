/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file omsetparm.cpp

  \brief Setting up parameters in optical model calculations

*/

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>

using namespace std;

#include "optical.h"
#include "physics.h"
#include "structur.h"

static void omSetCoulomb(int, double, double, double *, double *);

/**********************************************************/
/*      Setup Energy Dependent Parameters for OM Calc     */
/*      cms energy, wave number, Coulomb parameter of     */
/*      ejectile, excitation of residual nucleus          */
/**********************************************************/
void omSetEnergy(double e, int zz, double mu, CCdata *c) {
  c->energy = e;
  c->reduced_mass = mu;
  c->wavesq = c->energy * mu * 2.0 * AMUNIT / VLIGHTSQ / HBARSQ;

  if (zz == 0) {
    c->coulomb = c->coulomb_scat0 = 0.0;
  } else {
    omSetCoulomb(zz, mu, c->energy, &c->coulomb, &c->coulomb_scat0);
  }
}

/**********************************************************/
/*     Coulomb Parameters                                 */
/**********************************************************/
void omSetCoulomb(int zz, double mu, double e, double *yeta, double *sig0) {
  double y1, y2, y3, y4;

  *yeta = y1 = PERMITTIV * COULOMBSQ * zz * sqrt(AMUNIT * mu / (2.0 * e)) /
               VLIGHT / HBAR;
  y2 = y1 * y1;
  y3 = 16.0 + y2;
  y4 = y3 * y3;
  *sig0 = -y1 + y1 * log(y3) / 2. + 3.5 * atan(y1 / 4.) -
          (atan(y1) + atan(y1 / 2.) + atan(y1 / 3.)) -
          y1 *
              (1. + (y2 - 48.) / (30. * y4) +
               (y2 * y2 - 160. * y2 + 1280.) / (105. * y4 * y4)) /
              (12. * y3);
}

/**********************************************************/
/*     Setup Optical Potential Parameters                 */
/**********************************************************/
unsigned int omSetOmp(double e, ZAnumber *target, Pdata *proj, Optical *omp) {
  auto potfm = proj->omp;
  const auto z0 = target->getZ();
  const auto a0 = target->getA();
  const auto z1 = proj->particle.getZ();
  const auto a1 = proj->particle.getA();
  omInitOmp(omp);
  unsigned int potfm0 = potfm;

  double a3 = pow((double)a0, 1.0 / 3.0);

  if (auto omp_file = proj->omp_file) {
    
    //  set potfm
    potfm = omp_library(potfm >> 8, z0, a0, z1, a1, e, omp);
    
    omp->volume.real = omp_file->real_cent_V(z0,a0,e);
    omp->volume.imag = omp_file->cmpl_cent_V(z0, a0, e);
    omp->surface.real = omp_file->real_surf_V(z0, a0, e);
    omp->surface.imag = omp_file->cmpl_surf_V(z0, a0, e);
    omp->spin_orbit.real = omp_file->real_spin_V(z0, a0, e);
    omp->spin_orbit.imag = omp_file->cmpl_spin_V(z0, a0, e);

    // diffusivities
    omp->a0 = omp_file->real_cent_a(z0, a0, e);
    omp->av = omp_file->cmpl_cent_a(z0, a0, e);
    omp->a0s = omp_file->real_surf_a(z0, a0, e);
    omp->as = omp_file->cmpl_surf_a(z0, a0, e);
    omp->avso = omp_file->real_spin_a(z0, a0, e);
    omp->awso = omp_file->cmpl_spin_a(z0, a0, e);

    // regular potential radia 
    omp->R0 = omp_file->real_cent_r(z0, a0, e);
    omp->Rv = omp_file->cmpl_cent_r(z0, a0, e);
    omp->R0s = omp_file->real_surf_r(z0, a0, e);
    omp->Rs = omp_file->cmpl_surf_r(z0, a0, e);
    omp->Rvso = omp_file->real_spin_r(z0, a0, e);
    omp->Rwso = omp_file->cmpl_spin_r(z0, a0, e);
    omp->Rc = 0;
    
    // reduced radii
    omp->r0 = omp->R0 / a3;
    omp->r0s = omp->R0s /  a3;
    omp->rv = omp->Rv / a3;
    omp->rs = omp->Rs / a3;
    omp->rvso = omp->Rvso / a3;
    omp->rwso = omp->Rwso / a3;
    omp->rc = omp->Rc / a3;

  } else {
    //  Retrieve potential parameters from library
    potfm = omp_library(potfm >> 8, z0, a0, z1, a1, e, omp);

    //  Energy dependent depths
    omp->volume.real = omp->v1 + omp->v2 * e + omp->v3 * e * e;
    omp->surface.real = omp->vs1 + omp->vs2 * e + omp->vs3 * e * e;
    omp->volume.imag = omp->wv1 + omp->wv2 * e + omp->wv3 * e * e;
    omp->surface.imag = omp->ws1 + omp->ws2 * e + omp->ws3 * e * e;
    omp->spin_orbit.real = omp->vso1 + omp->vso2 * e + omp->vso3 * e * e;
    omp->spin_orbit.imag = omp->wso1 + omp->wso2 * e + omp->wso3 * e * e;
    
    // reduced potential radii
    omp->R0 = omp->r0 * a3;
    omp->R0s = omp->r0s * a3;
    omp->Rv = omp->rv * a3;
    omp->Rs = omp->rs * a3;
    omp->Rvso = omp->rvso * a3;
    omp->Rwso = omp->rwso * a3;
    omp->Rc = omp->rc * a3;
  }
  

  if (omp->volume.real <= 0.0) {
    potfm = (potfm & 0x009f);
    omp->volume.real = 0.0;
  }
  if (omp->volume.imag <= 0.0) {
    potfm = (potfm & 0x00fe);
    omp->volume.imag = 0.0;
  }
  if (omp->surface.imag <= 0.0) {
    potfm = (potfm & 0x00f9);
    omp->surface.imag = 0.0;
  }
  if (omp->spin_orbit.real <= 0.0) {
    omp->spin_orbit.real = 0.0;
  }

  return ((potfm0 & 0xff00) | potfm);
}

/**********************************************************/
/*     Clear Optical Potential Parameters                 */
/**********************************************************/
void omInitOmp(Optical *omp) {
  omp->r0 = omp->r0s = omp->rv = omp->rs = omp->rvso = omp->rwso = 0.0;
  omp->a0 = omp->a0s = omp->av = omp->as = omp->avso = omp->awso = 0.0;
  omp->v1 = omp->vs1 = omp->wv1 = omp->ws1 = omp->vso1 = omp->wso1 = 0.0;
  omp->v2 = omp->vs2 = omp->wv2 = omp->ws2 = omp->vso2 = omp->wso2 = 0.0;
  omp->v3 = omp->vs3 = omp->wv3 = omp->ws3 = omp->vso3 = omp->wso3 = 0.0;
  omp->rc = 0.0;
}
