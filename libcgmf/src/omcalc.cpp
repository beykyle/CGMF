/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file omcalc.cpp

  \brief Spherical optical model calculation

*/

#include <cmath>
#include <cstdio>
#include <iostream>
#include <cassert>

using namespace std;

#include "cgm.h"
#include "optical.h"
#include "structur.h"

/**********************************************************/
/*     Entrance channel calculation                       */
/*         usual optical model calculation for incident   */
/*         particle and determine max-J                   */
/**********************************************************/
int omCalc(double energy,    // CMS incident energy
           Pdata *proj,      // particle data for projectile
           ZAnumber *targ,   // target ZA
           double mu,        // reduced mass for Targ+Proj system
           double *tran,     // transmission coefficients (output)
           CrossSection *crx // calculated cross sections (output)
) {
  Optical omp;
  Potential pot;
  Wavefunc wfn;
  CCdata cdt;
  Complex *smat;

  int l, zzprod = 0;
  double cx[3];

  if (energy <= 0.0)
    return (-1);
  if (proj->particle.getZ() > 0 && energy < ECUT_CHARGED)
    return (-1);


  if (proj->emulate) {
    crx->elastic = crx->reaction = crx->total = 0.0;
    cdt.lmax = 10;
    Complex c, d;
    int jmax = (int)(2 * proj->spin) + 1;
    assert(proj->particleID == neutron);
    for (int j = 0; j < 3 * MAX_L; j++)

    zzprod = targ->getZ() * proj->particle.getZ();
    omSetEnergy(energy, zzprod, mu, &cdt);
    unsigned int ompindex = omSetOmp(cdt.energy, targ, proj, &omp) & 0x00ff;
    pot.width = INTEG_WIDTH;
    omPotentialForm(ompindex, zzprod, &omp, &cdt, &pot, true);
    
    // emulation


    for (l = 0; l <= cdt.lmax; l++) {
      Complex smat[jmax];
      for (int j = 0; j < jmax; j++) {
        double xj = l + proj->spin - (double)j;
        if (xj < fabs(l - proj->spin))
          continue;
        int index = l * 3 + j;
        smat[j] = Complex(0.146729937, 0.03935661271);
        tran[index] = 1.0 - absolute(&smat[j]);
        if (tran[index] < 0.0)
          tran[index] = 0.0;
      }
      
      c.real = 1 - smat[0].real;
      c.imag = smat[0].imag;
      d.real = 1 - smat[1].real;
      d.imag = smat[1].imag;
      crx->elastic += (l + 1) * absolute(&c) + l * absolute(&d);
      crx->reaction +=  (l + 1) * tran[l * 3] + l * tran[l * 3 + 1];
    }
    return cdt.lmax;
  }
  
  //---------------------------------------
  //      Memory Allocation
  try {
    pot.mean_field = new Complex[MAX_POINTS];
    pot.spin_orbit = new Complex[MAX_POINTS];
    pot.radi = new double[MAX_POINTS];
    pot.coulomb = new double[MAX_POINTS];
    wfn.external = new Complex[MAX_L];
    wfn.extderiv = new Complex[MAX_L];
    wfn.internal = new Complex[MAX_POINTS];
    smat = new Complex[MAX_L * 3];
  } catch (bad_alloc) {
    cerr << "Memory allocation error";
    return (-1);
  }

  for (int j = 0; j < 3 * MAX_L; j++)
    tran[j] = 0.0;

  zzprod = targ->getZ() * proj->particle.getZ();

  //---------------------------------------
  //      Energy and Coulomb Parameter

  omSetEnergy(energy, zzprod, mu, &cdt);

  //---------------------------------------
  //      Setup Optical Potential Geometry

  unsigned int ompindex = omSetOmp(cdt.energy, targ, proj, &omp) & 0x00ff;
  pot.width = INTEG_WIDTH;
  omPotentialForm(ompindex, zzprod, &omp, &cdt, &pot);

  //---------------------------------------
  //      Main Calculation

  /***  free space wave function */
  cdt.lmax = omExternalFunction(0, pot.rho_match, cdt.coulomb,
                                cdt.coulomb_scat0, &wfn);

  /***  internal wave function */
  Complex a, b, c, d;
  int jmax = (int)(2 * proj->spin) + 1;

  crx->elastic = crx->reaction = crx->total = 0.0;
  for (l = 0; l <= cdt.lmax; l++) {
    a = rational(wfn.extderiv[l].real, wfn.extderiv[l].imag,
                 wfn.external[l].real, wfn.external[l].imag);
    b = rational(wfn.external[l].real, -wfn.external[l].imag,
                 wfn.external[l].real, wfn.external[l].imag);

    a.real *= pot.rho_match;
    a.imag *= pot.rho_match;

    for (int j = 0; j < jmax; j++) {
      double xj = l + proj->spin - (double)j;
      
      if (xj < fabs(l - proj->spin))
        continue;
      omInternalFunction(pot.n_match, pot.width, cdt.wavesq, (double)l,
                         proj->spin, xj, &pot, &wfn);

      int index = l * 3 + j;
      smat[index] =
          omSmatrix(pot.n_match, pot.width, pot.rad_match, &a, &b, &wfn);
      tran[index] = 1.0 - absolute(&smat[index]);
      if (tran[index] < 0.0)
        tran[index] = 0.0;
    }

    /***  reaction cross section */
    cx[0] = crx->reaction;
    cx[1] = 0.0;
    switch (proj->particleID) {
    case neutron:
      c.real = 1 - smat[l * 3].real;
      c.imag = smat[l * 3].imag;
      d.real = 1 - smat[l * 3 + 1].real;
      d.imag = smat[l * 3 + 1].imag;
      crx->elastic += (l + 1) * absolute(&c) + l * absolute(&d);
    case proton:
    case triton:
    case helion:
      crx->reaction += (cx[1] = (l + 1) * tran[l * 3] + l * tran[l * 3 + 1]);
      break;
    case alpha:
      crx->reaction += (cx[1] = (2 * l + 1) * tran[l * 3]);
      break;
    case deuteron:
      crx->reaction +=
          (cx[1] = ((2 * l + 3) * tran[l * 3] + (2 * l + 1) * tran[l * 3 + 1] +
                    (2 * l - 1) * tran[l * 3 + 2]) /
                   3.0);
      break;
    default:
      break;
    }

    if (cx[0] != 0.0 && (cx[1] / cx[0]) < CRIT_LCUT)
      break;
  }
  if (l == MAX_L - 1) {
    cerr << "reaction cross section not converge" << endl;
    cdt.lmax = l;
  } else if (l < cdt.lmax)
    cdt.lmax = l;

  //---------------------------------------
  //      Output Cross Sections

  //  crx->elastic  *= NORM_FACT*PI/cdt.wavesq;
  //  crx->reaction *= NORM_FACT*PI/cdt.wavesq;
  //  crx->total  = crx->elastic + crx->reaction;

  //---------------------------------------
  //      Free Allocated Memory

  delete[] pot.mean_field;
  delete[] pot.spin_orbit;
  delete[] pot.radi;
  delete[] pot.coulomb;
  delete[] wfn.external;
  delete[] wfn.extderiv;
  delete[] wfn.internal;
  delete[] smat;

  return (cdt.lmax);
}
