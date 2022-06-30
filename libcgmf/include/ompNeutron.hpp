#ifndef __GKDN__
#define __GKDN__

#ifdef __cplusplus
#include <string>
using std::string;
#include "json.hpp"
using nlohmann::json;
#endif

struct OMPFile {
  virtual double real_central_radius(int zt, int at, double En) const = 0;
  virtual double cmpl_central_radius(int zt, int at, double En) const = 0;
  virtual double real_surf_radius(int zt, int at, double En) const = 0;
  virtual double cmpl_surf_radius(int zt, int at, double En) const = 0;
  virtual double real_so_radius(int zt, int at, double En) const   = 0;
  virtual double cmpl_so_radius(int zt, int at, double En) const   = 0;
  
  virtual double real_central_diffusivity(int zt, int at, double En) const = 0;
  virtual double cmpl_central_diffusivity(int zt, int at, double En) const = 0;
  virtual double real_surf_diffusivity(int zt, int at, double en) const = 0;
  virtual double cmpl_surf_diffusivity(int zt, int at, double En) const = 0;
  virtual double real_so_diffusivity(int zt, int at, double En) const = 0;
  virtual double cmpl_so_diffusivity(int zt, int at, double En) const = 0;

  virtual double real_central_depth(int zt, int at, double En) const = 0;
  virtual double cmpl_central_depth(int zt, int at, double En) const = 0;
  virtual double real_surf_depth(int zt, int at, double En) const  = 0;
  virtual double cmpl_surf_depth(int zt, int at, double En) const = 0;
  virtual double real_so_depth(int zt, int at, double En) const = 0;
  virtual double cmpl_so_depth(int zt, int at, double En) const = 0;
};

struct GKDNeutron : public OMPFile {
  const double e_fermi_0;
  const double e_fermi_A;
  
  // real central
  const double v1_0, v1_asym, v1_A, v2_0, v2_A, v3_0, v3_A, v4_0; // depth
  const double r_0, r_A, a_0, a_A;                                // shape

  // cmplex central
  const double w1_0, w1_A, w2_0, w2_A;
  
  // cmplex surface
  const double d1_0, d1_asym, d2_0, d2_A, d2_A2, d2_A3, d3_0;
  const double rd_0, rd_A, ad_0, ad_A;

  // real spin orbit
  const double vso1_0, vso1_A, vso2_0;
  const double rso_0, rso_A, aso_0;
  
  // cmplex spin orbit
  const double wso1, wso2;

  // structure and energy factors
  double Ef(int at) const;
  double asym(int zt, int at) const;

  // pot radii
  double real_central_radius(int zt, int at, double En) const final;
  double cmpl_central_radius(int zt, int at, double En) const final;
  double real_surf_radius(int zt, int at, double En) const final { return 0; };
  double cmpl_surf_radius(int zt, int at, double En) const final;
  double real_so_radius(int zt, int at, double En) const final;
  double cmpl_so_radius(int zt, int at, double En) const final;
  
  // pot diffusivities
  double real_central_diffusivity(int zt, int at, double En) const final;
  double cmpl_central_diffusivity(int zt, int at, double En) const final;
  double real_surf_diffusivity(int zt, int at, double en) const final { return 0; }
  double cmpl_surf_diffusivity(int zt, int at, double En) const final;
  double real_so_diffusivity(int zt, int at, double En) const final;
  double cmpl_so_diffusivity(int zt, int at, double En) const final;

  // pot depths
  double real_central_depth(int zt, int at, double En) const final;
  double cmpl_central_depth(int zt, int at, double En) const final;
  double real_surf_depth(int zt, int at, double En) const final {return 0;};
  double cmpl_surf_depth(int zt, int at, double En) const final;
  double real_so_depth(int zt, int at, double En) const final;
  double cmpl_so_depth(int zt, int at, double En) const final; 
  
  // read params from json file
  GKDNeutron(json param_file);

  // set default KD global params
  GKDNeutron(): 
    // real central
    v1_0(59.30)   , v1_asym(21.0), v1_A(0.024)    , v2_0( 0.007228), v2_A(1.48e-6)
  , v3_0(1.994e-5), v3_A( 2.0e-8), v4_0(7e-9)     , r_0(1.3039)    , r_A(0.4054)    
  , a_0(0.6778)   , a_A(1.487e-4) 
  // cmplex central
  , w1_0(12.195)  , w1_A(0.0167) , w2_0(73.55)    , w2_A(0.0795)
  // cmplex surface 
  , d1_0(16.0)    , d1_asym(16.0), d2_0(0.0180)   , d2_A(0.003802) , d2_A2(8.0)
  , d2_A3(156.0)  , d3_0(11.5)   , rd_0(1.3424)   , rd_A(0.01585)  , ad_0(0.5446)
  , ad_A(1.656e-4)
  // real spin orbit
  , vso1_0(5.922) , vso1_A(0.0030), vso2_0(0.0040), rso_0(1.1854)  , rso_A(0.647)
  , aso_0(0.59)
  // cmplex spin orbit
  , wso1(-3.1)    , wso2(160) 
  // fermi energy
  , e_fermi_0(-11.2814), e_fermi_A(0.02646)
  {}
};

struct WLHNeutron : public OMPFile {
  // pot radii
  double real_central_radius(int zt, int at, double En) const final;
  double cmpl_central_radius(int zt, int at, double En) const final;
  double real_surf_radius(int zt, int at, double En) const final { return 0; };
  double cmpl_surf_radius(int zt, int at, double En) const final;
  double real_so_radius(int zt, int at, double En) const final;
  double cmpl_so_radius(int zt, int at, double En) const final;
  
  // pot diffusivities
  double real_central_diffusivity(int zt, int at, double En) const final;
  double cmpl_central_diffusivity(int zt, int at, double En) const final;
  double real_surf_diffusivity(int zt, int at, double en) const final { return 0; }
  double cmpl_surf_diffusivity(int zt, int at, double En) const final;
  double real_so_diffusivity(int zt, int at, double En) const final;
  double cmpl_so_diffusivity(int zt, int at, double En) const final;

  // pot depths
  double real_central_depth(int zt, int at, double En) const final;
  double cmpl_central_depth(int zt, int at, double En) const final;
  double real_surf_depth(int zt, int at, double En) const final {return 0;};
  double cmpl_surf_depth(int zt, int at, double En) const final;
  double real_so_depth(int zt, int at, double En) const final;
  double cmpl_so_depth(int zt, int at, double En) const final; 

  WLHNeutron(json param_file);
};

struct CH89Neutron : public OMPFile {
  // pot radii
  double real_central_radius(int zt, int at, double En) const final;
  double cmpl_central_radius(int zt, int at, double En) const final;
  double real_surf_radius(int zt, int at, double En) const final { return 0; };
  double cmpl_surf_radius(int zt, int at, double En) const final;
  double real_so_radius(int zt, int at, double En) const final;
  double cmpl_so_radius(int zt, int at, double En) const final;
  
  // pot diffusivities
  double real_central_diffusivity(int zt, int at, double En) const final;
  double cmpl_central_diffusivity(int zt, int at, double En) const final;
  double real_surf_diffusivity(int zt, int at, double en) const final { return 0; }
  double cmpl_surf_diffusivity(int zt, int at, double En) const final;
  double real_so_diffusivity(int zt, int at, double En) const final;
  double cmpl_so_diffusivity(int zt, int at, double En) const final;

  // pot depths
  double real_central_depth(int zt, int at, double En) const final;
  double cmpl_central_depth(int zt, int at, double En) const final;
  double real_surf_depth(int zt, int at, double En) const final {return 0;};
  double cmpl_surf_depth(int zt, int at, double En) const final;
  double real_so_depth(int zt, int at, double En) const final;
  double cmpl_so_depth(int zt, int at, double En) const final; 
  
  CH89Neutron(json param_file);
};

#endif
