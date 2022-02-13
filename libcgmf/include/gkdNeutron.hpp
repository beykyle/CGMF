#ifndef __GKDN__
#define __GKDN__

#ifdef __cplusplus
#include <string>
using std::string;
#endif

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

  // structure factor 
  double asym(int zt, int at) const;

  // potential shapes
  double real_radius(int zt, int at, double e) const;
  double so_radius(int zt, int at, double e) const;
  double compl_surf_radius(int zt, int at, double e) const;
  double real_diffusivity(int zt, int at, double e) const;
  double so_diffusivity(int zt, int at, double e) const;
  double compl_surf_diffusivity(int zt, int at, double e) const;

  // potential depths
  double real_central_depth(int zt, int at, double e) const;
  double compl_central_depth(int zt, int at, double e) const;
  double compl_surf_depth(int zt, int at, double e) const;
  double real_so_depth(int zt, int at, double e) const;
  double compl_so_depth(int zt, int at, double e) const; 

  // read params from json file
  GKDNeutron(string fname);

  // set default KD global params
  GKDNeutron(): 
    // real central
    v1_0(59.30)   , v1_asym(21.0), v1_A(0.024)    , v2_0( 0.007228), v2_A(1.48e-6)
  , v3_0(1.994e-5), v3_A( 2.0e-8), v4_0(7e-9)     , r_0(1.3039)    , r_A(0.4054)    
  , a_0(0.6778)   , a_A(1.487e-4) 
  // complex central
  , w1_0(12.195)  , w1_A(0.0167) , w2_0(73.55)    , w2_A(0.0795)
  // complex surface 
  , d1_0(16.0)    , d1_asym(16.0), d2_0(0.0180)   , d2_A(0.003802) , d2_A2(8.0)
  , d2_A3(156.0)  , d3_0(11.5)   , rd_0(1.3424)   , rd_A(0.01585)  , ad_0(0.05446)
  , ad_A(1.656e-4)
  // real spin orbit
  , vso1_0(5.922) , vso1_A(0.0030), vso2_0(0.0040), rso_0(1.1854)  , rso_A(0.647)
  , aso_0(0.59)
  // complex spin orbit
  , wso1(-3.1)    , wso2(160) {}
};

#endif
