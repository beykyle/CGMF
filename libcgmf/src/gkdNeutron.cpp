#include <iostream>
#include <cmath>

using namespace std;

#include "physics.h"
#include "gkdNeutron.hpp"
#include "terminate.h"


double GKDNeutron::asym(int zt, int at) const {
  const double A = (double)at;
  const double Z = (double)zt;
  return 1 - 2*A/Z;
}

double GKDNeutron::real_radius(int zt, int at, double e) const {
  const double A = (double)at;
  return r_0 - r_A * pow(A,-1./3.);
}

double GKDNeutron::so_radius(int zt, int at, double e) const {
  const double A = (double)at;
  return rso_0 - rso_A * pow(A,-1./3.);
}

double GKDNeutron::compl_surf_radius(int zt, int at, double e) const {
  const double A = (double)at;
  return rd_0 - rd_A * pow(A,1./3.);
}


double GKDNeutron::real_diffusivity(int zt, int at, double e) const {
  const double A = (double)at;
  return a_0 - a_A * A;
}

double GKDNeutron::so_diffusivity(int zt, int at, double e) const {
  return aso_0;
}

double GKDNeutron::compl_surf_diffusivity(int zt, int at, double e) const {
  const double A = (double)at;
  return ad_0 - ad_A * A;
}

double GKDNeutron::real_central_depth(int zt, int at, double e) const {
  const double Ex = e - e_fermi;
  const double A = (double)at;
  const double alpha = asym(zt,at);
  
  const double v1 = v1_0 - v1_asym * alpha - v1_A * A;
  const double v2 = v2_0 - v2_A * A;
  const double v3 = v3_0 - v3_A * A;
  const double v4 = v4_0;

  return v1 * (1 -  v2 * Ex + v3 * Ex*Ex - v4 * Ex*Ex*Ex);
}

double GKDNeutron::compl_central_depth(int zt, int at, double e) const {
  const double Ex = e - e_fermi;
  const double A = (double)at;

  const double w1 = w1_0 + w1_A * A;
  const double w2 = w2_0 + w2_A * A;

  return w1 * Ex * Ex/(Ex*Ex + w2*w2);
}

double GKDNeutron::compl_surf_depth(int zt, int at, double e) const {
  const double Ex = e - e_fermi;
  const double A = (double)at;
  const double alpha = asym(zt,at);

  const double d1 = d1_0 - d1_asym * alpha;
  const double d2 = d2_0 + d2_A /(1 +  exp( (A - d2_A3)/d2_A2) );
  const double d3 = d3_0;

  return d1 * Ex * Ex/(Ex*Ex + d3*d3) * exp( -d2 * Ex);
}

double GKDNeutron::real_so_depth(int zt, int at, double e) const {
  const double Ex = e - e_fermi;
  const double A = (double)at;

  const double vso1 = vso1_0 + vso1_A * A;
  const double vso2 = vso2_0;

  return vso1 * exp( -vso2 * Ex);
}

double GKDNeutron::compl_so_depth(int zt, int at, double e) const {
  const double Ex = e - e_fermi;
  return wso1 * Ex * Ex/(Ex*Ex + wso2*wso2);
}

GKDNeutron::GKDNeutron(string fname) {
  //TODO
}
