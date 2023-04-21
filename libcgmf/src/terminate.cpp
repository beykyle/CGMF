/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file terminate.cpp

  \brief Stops & Exits

*/

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;

#include "cgm.h"
#include "terminate.h"

/**********************************************************/
/*     Emergency Stop                                     */
/**********************************************************/
int cgmTerminateCode(std::string msg) {
  /*** Release global storage */
  cgmDeleteAllocated();

  /*** Exit code */
  cerr << "ERROR     :" << msg << endl;
  exit(-1);
}

int cgmTerminateCode(std::string msg, int n) {
  /*** Release global storage */
  cgmDeleteAllocated();

  /*** Exit code */
  cerr << "ERROR     :" << msg << " : " << n << endl;
  exit(-1);
}

int cgmTerminateCode(std::string msg, double x) {
  /*** Release global storage */
  cgmDeleteAllocated();

  /*** Exit code */
  cerr << "ERROR     :" << msg << " : " << x << endl;
  exit(-1);
}

/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void cgmDeleteAllocated() {
  delete[] ncl;
  for (int i = 0; i < SPECTRA_OUTPUT; i++)
    delete[] spc[i];
}

/**********************************************************/
/*     Memory Allocation                                  */
/**********************************************************/
void cgmAllocateMemory() {
  try {
    /*** compound nucleus */
    ncl = new Nucleus[MAX_COMPOUND];

    /*** calculated results */
    for (int i = 0; i < SPECTRA_OUTPUT; i++) {
      spc[i] = new double[MAX_ENERGY_BIN];
    }
  } catch (bad_alloc) {
    cerr << "ERROR     :memory allocation error";
    exit(-1);
  }
}
