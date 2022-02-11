/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file optical.h

  \brief Global definitions, constants and function prototypes for optical model calculations

*/
#ifndef __OPTICAL_H__
#define __OPTICAL_H__

#include "maths.h"

/**************************************/
/*      Constants                     */
/**************************************/
#define MAX_ITER        100         /* maximum iteration number            */
#define MAX_L            50         /* maximum angular momentum + 2        */
#define MAX_L0          100         /* maximum iteration for F function    */
#define MAX_LTRANS       10         /* maximum angular momentum transfer   */

#define MAX_POINTS      300         /* maximum radial points for unbound wave */
#define INTEG_WIDTH     0.2         /* default wave func. integral width   */
#define CRIT_MATCHING   1.0e-06     /* matching radius criterion           */
#define CRIT_LMAX       1.0e-10     /* maximum angular momentum criterion  */
#define CRIT_LCUT       1.0e-06     /* reaction cross section converge     */
#define CRIT_SCUT       1.0e-06     /* spin distribution cut off           */
#define CRIT_WRONSK     1.0e-06     /* wronskian satisfaction criterion    */
#define CRIT_ITER       1.0e-16     /* infinit iteration criterion         */
#define CRIT_LEVDEN     1.0e-06     /* level density matching criterion    */

#define ECUT_CHARGED    1.0e-10     /* cut-off energy for chared particles */


/**************************************/
/*      Structure data                */
/**************************************/

/****************************/
/*   OM Potential Parameter */
/****************************/
struct GKDNeutron {
  double e_fermi;
  // real central
  double v1_0, v1_asym, v1_A, v2_0, v2_A, v3_0, v3_A, v4_0; // depth
  double r_0, r_A, a_0, a_A;                                // shape

  // complex central
  double w1_0, w1_A, w2_0, w2_A;
  
  // complex surface
  double d1_0, d1_asym, d2_0, d2_A, d2_A2, d2_A3, d3_0;
  double rd_0, rd_A, ad_0, ad_A;

  // real spin orbit
  double vso1_0, vso1_A, vso2_0;
  double rso_0, rso_A, aso_0;
  
  // complex spin orbit
  double wso1, wso2;

  double asym(int zt, int at) const {
    const double A = (double)at;
    const double Z = (double)zt;
    return 1 - 2*A/Z;
  }

  double real_radius(int zt, int at, double e) const {
    const double A = (double)at;
    return r_0 - r_A * pow(A,1./3.);
  }
  
  double so_radius(int zt, int at, double e) const {
    const double A = (double)at;
    return rso_0 - rso_A * pow(A,1./3.);
  }
  
  double compl_surf_radius(int zt, int at, double e) const {
    const double A = (double)at;
    return rd_0 - rd_A * pow(A,1./3.);
  }
  
  double real_diffusivity(int zt, int at, double e) const {
    const double A = (double)at;
    return a_0 - a_A * A;
  }
  
  double so_diffusivity(int zt, int at, double e) const {
    return aso_0;
  }
  
  double compl_surf_diffusivity(int zt, int at, double e) const {
    const double A = (double)at;
    return ad_0 - ad_A * A;
  }

  double real_central_depth(int zt, int at, double e) const {
    const double Ex = e - e_fermi;
    const double A = (double)at;
    const double alpha = asym(zt,at);
    
    const double v1 = v1_0 - v1_asym * alpha - v1_A * A;
    const double v2 = v2_0 - v2_A * A;
    const double v3 = v3_0 - v3_A * A;
    const double v4 = v4_0;

    return v1 * (1 - 2 * v2 * Ex + v3 * Ex*Ex - v4 * Ex*Ex*Ex);
  }

  double compl_central_depth(int zt, int at, double e) const {
    const double Ex = e - e_fermi;
    const double A = (double)at;

    const double w1 = w1_0 + w1_A * A;
    const double w2 = w2_0 + w2_A * A;

    return w1 * Ex * Ex/(Ex*Ex + w2*w2);
  }
  
  double compl_surf_depth(int zt, int at, double e) const {
    const double Ex = e - e_fermi;
    const double A = (double)at;
    const double alpha = asym(zt,at);

    const double d1 = d1_0 - d1_asym * alpha;
    const double d2 = d2_0 + d2_A /(1 +  exp( (A - d2_A3)/d2_A2) );
    const double d3 = d3_0;

    return d1 * Ex * Ex/(Ex*Ex + d3*d3) * exp( -d2 * Ex);
  }

  double real_so_depth(int zt, int at, double e) const {
    const double Ex = e - e_fermi;
    const double A = (double)at;

    const double vso1 = vso1_0 + vso1_A * A;
    const double vso2 = vso2_0;

    return vso1 * exp( -vso2 * Ex);
  }
  
  double compl_so_depth(int zt, int at, double e) const {
    const double Ex = e - e_fermi;
    return wso1 * Ex * Ex/(Ex*Ex + wso2*wso2);
  }


  GKDNeutron(int zt, int at, string fname) {};
};

class Optical{
 public:
    double r0  ;double r0s ;double rv  ;double rs  ;double rvso;double rwso;
    double a0  ;double a0s ;double av  ;double as  ;double avso;double awso;
    double v1  ;double vs1 ;double wv1 ;double ws1 ;double vso1;double wso1;
    double v2  ;double vs2 ;double wv2 ;double ws2 ;double vso2;double wso2;
    double v3  ;double vs3 ;double wv3 ;double ws3 ;double vso3;double wso3;
    double rc  ;
    double R0  ;double R0s ;double Rv  ;double Rs  ;double Rvso;double Rwso;
    double Rc  ;
    Complex volume    ;
    Complex surface   ;
    Complex spin_orbit;

    Optical(){
     r0  = 0.0; r0s = 0.0; rv  = 0.0; rs  = 0.0; rvso= 0.0; rwso= 0.0;
     a0  = 0.0; a0s = 0.0; av  = 0.0; as  = 0.0; avso= 0.0; awso= 0.0;
     v1  = 0.0; vs1 = 0.0; wv1 = 0.0; ws1 = 0.0; vso1= 0.0; wso1= 0.0;
     v2  = 0.0; vs2 = 0.0; wv2 = 0.0; ws2 = 0.0; vso2= 0.0; wso2= 0.0;
     v3  = 0.0; vs3 = 0.0; wv3 = 0.0; ws3 = 0.0; vso3= 0.0; wso3= 0.0;
     rc  = 0.0;
     R0  = 0.0; R0s = 0.0; Rv  = 0.0; Rs  = 0.0; Rvso= 0.0; Rwso= 0.0;
     Rc  = 0.0;
    }
};

/****************************/
/*   Optical Potential      */
/****************************/
class Potential{
 public:
    int       n_match        ;     /* integration step number          */
    double    rad_match      ;     /* matching radius                  */
    double    rho_match      ;     /* rad_match * wave_number          */
    double    width          ;     /* integration step width           */
    double  *radi            ;     /* r*r                              */
    double  *coulomb         ;     /* Coulomb field                    */
    Complex  *mean_field     ;     /* Optical potential V   + iW       */
    Complex *spin_orbit      ;     /* Optical potenital Vso + iWso     */

    Potential(){
      n_match   = 0;
      rad_match = 0.0;
      rho_match = 0.0;
      width     = 0.0;
    }
};

/****************************/
/*   Unbound State Data     */
/****************************/
class Wavefunc{
 public:
    Complex *external        ;     /* external wave func.  G, F at Rm  */
    Complex *extderiv        ;     /* derivative ext.wave dG,dF at Rm  */
    Complex *internal        ;     /* internal wave func. R=0 to Rm    */
};


/****************************/
/*   Channel Spin           */
/****************************/
class Chspin{
 public:
    int      l               ;     /* angular momentum for channel     */
    int      j2              ;     /* 2 x j (j = l+s )                 */
    int      pt              ;     /* parity of target                 */
    int      pi              ;     /* parity of emitted particle       */
    double   si              ;     /* spin of emitted particle         */ 
    double   st              ;     /* spin of target                   */

    Chspin(){
      l  = 0;
      j2 = 0;
      pt = 0;
      pi = 0;
      si = 0.0;
      st = 0.0;
    }
};

/****************************/
/*   Channel Data           */
/****************************/
class CCdata{
 public:
    Chspin   chn             ;     /* channel spin                     */
    int      lmax            ;     /* max angular momentum for channel */
    int      jmax            ;     /* max spin                         */
    int      level           ;     /* index label for excited states   */
    double   energy          ;     /* energy of the emitte particle    */
    double   wavesq          ;     /* square of wave number            */
    double   reduced_mass    ;     /* reduced mass                     */
    double   coulomb         ;     /* Coulomb parameter                */
    double   coulomb_scat0   ;     /* Coulomb phase shift for L=0      */
    double   xl_term         ;     /* L(L+1)                           */
    double   so_term         ;     /* J(J+1)-L(L+1)-S(S+1)             */

    CCdata(){
      lmax          = 0;
      jmax          = 0;
      level         = 0;
      energy        = 0.0;
      wavesq        = 0.0;
      reduced_mass  = 0.0;
      coulomb       = 0.0;
      coulomb_scat0 = 0.0;
      xl_term       = 0.0;
      so_term       = 0.0;
    }
};


/**************************************/
/*      asympt.cpp                    */
/**************************************/
void    omAsymptotic           (double, double, double, double *, double *);


/**************************************/
/*      builtin.cpp                   */
/**************************************/
unsigned int  omp_library      (int, int, int, int, int, double, Optical *);
unsigned int  find_omp         (string);


/**************************************/
/*      extwave.cpp                   */
/**************************************/
int     omExternalFunction     (int, double, double, double, Wavefunc *);
int     omExternalClosed       (int, double, double, Wavefunc *);


/**************************************/
/*      intwave.cpp                   */
/**************************************/
void    omInternalFunction     (int, double, double, double, double, double,
                                Potential *, Wavefunc *);
void    omIrregularFunction    (int, double, double, double, double, double,
                                Potential *, Wavefunc *);
Complex omNormalizationFactor  (int, int, double, double, Complex *,Wavefunc *);
Complex omCoulombPhaseFactor   (int,      double, double);
void    omNormalization        (int,      double, Complex *,Potential *,Wavefunc *);


/**************************************/
/*      omsetform.cpp                 */
/**************************************/
int     omPotentialForm        (unsigned int, int, Optical *, CCdata *,Potential *);
int     omPotentialFixedLength (int, unsigned int, int, Optical *, CCdata *,Potential *);
double  omPotentialRadialCoulomb(double, Optical *);


/**************************************/
/*      omsetparm.cpp                 */
/**************************************/
void    omSetEnergy            (double, int, double, CCdata *);
void    omInitOmp              (Optical *);
unsigned int omSetOmp          (unsigned int, double, int, int, int, int, Optical *);


/**************************************/
/*      smatrix.cpp                   */
/**************************************/
Complex omSmatrix              (int, double,double,Complex *,Complex *,Wavefunc *);

#endif //__OPTICAL_H__
