#ifndef PYBINDINGS
#define PYBINDINGS

#include <string>

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "cgm.h"
#include "cgmfEvents.h"
#include "rngcgm.h"

#include "config-ff.h"
#include "config.h"

#ifdef MPIRUN
#include <mpi.h>
#endif

std::string test_func() {
  return "pyCGMF test";
}


namespace py = pybind11;
using arr = py::array_t<double>;

void run_histories_get_nubar( arr& nu) {
  assert(  nu.request().ndim != 1  );
  int nevents = nu.size();
  auto n = nu.mutable_unchecked<1>();
  const auto ZAIDt = 98252;
  const auto incidentEnergy = 0.0;
  const auto seed = 13;
  UniformRNG rng(1);
  rng.set_seed(seed);
  cgmfEvent* event = NULL;

  for ( int i = 0; i < nevents; i++ ) {
    event = new cgmfEvent(ZAIDt, incidentEnergy, 0.0, 1.0e-8, -1);
    auto nuLF = event->getLightFragmentNeutronNu();
    auto nuHF = event->getHeavyFragmentNeutronNu();
    auto nuPre = event->getPreFissionNeutronNu();
    n[i] = nuLF + nuHF + nuPre;
  }
}

//namespace py = pybind11;

PYBIND11_MODULE(pyCGMF, m) {
  m.doc() = "pyCGMF: python bindings for running CGMF on MPI";
  m.def(
      "run_histories_get_nubar", 
      &run_histories_get_nubar, 
      "runs nevents CGMF histories of 252-Cf (sf), returns the neutron multiplicity for each event" );
}

#endif
