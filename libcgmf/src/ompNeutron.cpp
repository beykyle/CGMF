#include <cmath>
#include <iostream>

using namespace std;

#include "ompNeutron.hpp"
#include "physics.h"
#include "terminate.h"

double GKDNeutron::Ef(int at) const {
  const double A = (double)at;
  return e_fermi_0 + e_fermi_A * A;
}

double GKDNeutron::asym(int zt, int at) const {
  const double A = (double)at;
  const double Z = (double)zt;
  return 1 - 2 * Z / A;
}

double GKDNeutron::real_central_radius(int zt, int at, double En) const {
  const double A = (double)at;
  return r_0 - r_A / pow(A, 1. / 3.);
}

double GKDNeutron::cmpl_central_radius(int zt, int at, double En) const {
  return real_central_radius(zt, at, En);
}

double GKDNeutron::real_so_radius(int zt, int at, double En) const {
  const double A = (double)at;
  return rso_0 - rso_A * pow(A, -1. / 3.);
}

double GKDNeutron::cmpl_so_radius(int zt, int at, double En) const {
  return real_so_radius(zt, at, En);
}

double GKDNeutron::cmpl_surf_radius(int zt, int at, double En) const {
  const double A = (double)at;
  return rd_0 - rd_A * pow(A, 1. / 3.);
}

double GKDNeutron::real_central_diffusivity(int zt, int at, double En) const {
  const double A = (double)at;
  return a_0 - a_A * A;
}

double GKDNeutron::cmpl_central_diffusivity(int zt, int at, double En) const {
  return real_central_diffusivity(zt, at, En);
}

double GKDNeutron::real_so_diffusivity(int zt, int at, double En) const {
  return aso_0;
}

double GKDNeutron::cmpl_so_diffusivity(int zt, int at, double En) const {
  return real_so_diffusivity(zt, at, En);
}

double GKDNeutron::cmpl_surf_diffusivity(int zt, int at, double En) const {
  const double A = (double)at;
  return ad_0 - ad_A * A;
}

double GKDNeutron::real_central_depth(int zt, int at, double En) const {
  const double Ex = En - Ef(at);
  const double A = (double)at;
  const double alpha = asym(zt, at);

  const double v1 = v1_0 - v1_asym * alpha - v1_A * A;
  const double v2 = v2_0 - v2_A * A;
  const double v3 = v3_0 - v3_A * A;
  const double v4 = v4_0;

  return v1 * (1 - v2 * Ex + v3 * Ex * Ex - v4 * Ex * Ex * Ex);
}

double GKDNeutron::cmpl_central_depth(int zt, int at, double En) const {
  const double Ex = En - Ef(at);
  const double A = (double)at;

  const double w1 = w1_0 + w1_A * A;
  const double w2 = w2_0 + w2_A * A;

  return w1 * Ex * Ex / (Ex * Ex + w2 * w2);
}

double GKDNeutron::cmpl_surf_depth(int zt, int at, double En) const {
  const double Ex = En - Ef(at);
  const double A = (double)at;
  const double alpha = asym(zt, at);

  const double d1 = d1_0 - d1_asym * alpha;
  const double d2 = d2_0 + d2_A / (1 + exp((A - d2_A3) / d2_A2));
  const double d3 = d3_0;

  return d1 * Ex * Ex / (Ex * Ex + d3 * d3) * exp(-d2 * Ex);
}

double GKDNeutron::real_so_depth(int zt, int at, double En) const {
  const double Ex = En - Ef(at);
  const double A = (double)at;

  const double vso1 = vso1_0 + vso1_A * A;
  const double vso2 = vso2_0;

  return vso1 * exp(-vso2 * Ex);
}

double GKDNeutron::cmpl_so_depth(int zt, int at, double En) const {
  const double Ex = En - Ef(at);
  return wso1 * Ex * Ex / (Ex * Ex + wso2 * wso2);
}

GKDNeutron::GKDNeutron(json p)
    : // real central
      v1_0(p["KDHartreeFock"]["V1_0"]), v1_asym(p["KDHartreeFock"]["V1_asymm"]),
      v1_A(p["KDHartreeFock"]["V1_A"]), v2_0(p["KDHartreeFock"]["V2_0_n"]),
      v2_A(p["KDHartreeFock"]["V2_A_n"]), v3_0(p["KDHartreeFock"]["V3_0_n"]),
      v3_A(p["KDHartreeFock"]["V3_A_n"]), v4_0(p["KDHartreeFock"]["V4_0"]),
      r_0(p["KDHartreeFock"]["r_0"]), r_A(p["KDHartreeFock"]["r_A"]),
      a_0(p["KDHartreeFock"]["a_0"]), a_A(p["KDHartreeFock"]["a_A"]),
      w1_0(p["KDImagVolume"]["W1_0_n"]), w1_A(p["KDImagVolume"]["W1_A_n"]),
      w2_0(p["KDImagVolume"]["W2_0"]), w2_A(p["KDImagVolume"]["W2_A"]),
      d1_0(p["KDImagSurface"]["D1_0"]), d1_asym(p["KDImagSurface"]["D1_asymm"]),
      d2_0(p["KDImagSurface"]["D2_0"]), d2_A(p["KDImagSurface"]["D2_A"]),
      d2_A2(p["KDImagSurface"]["D2_A2"]), d2_A3(p["KDImagSurface"]["D2_A3"]),
      d3_0(p["KDImagSurface"]["D3_0"]), rd_0(p["KDImagSurface"]["r_0"]),
      rd_A(p["KDImagSurface"]["r_A"]), ad_0(p["KDImagSurface"]["a_0_n"]),
      ad_A(p["KDImagSurface"]["a_A_n"]), vso1_0(p["KDRealSpinOrbit"]["V1_0"]),
      vso1_A(p["KDRealSpinOrbit"]["V1_A"]),
      vso2_0(p["KDRealSpinOrbit"]["V2_0"]), rso_0(p["KDRealSpinOrbit"]["r_0"]),
      rso_A(p["KDRealSpinOrbit"]["r_A"]), aso_0(p["KDRealSpinOrbit"]["a_0"]),
      wso1(p["KDImagSpinOrbit"]["W1_0"]), wso2(p["KDImagSpinOrbit"]["W2_0"]),
      e_fermi_0(-11.2814), e_fermi_A(0.02646) {}
