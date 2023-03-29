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

constexpr static auto ZAIDt = 98252;
constexpr static auto incidentEnergy = 0.0;
constexpr static auto seed = 13;

void run_histories_get_nubar( arr& nu) {
  // prepare array
  assert(  nu.request().ndim != 1  );
  int nevents = nu.size();
  auto n = nu.mutable_unchecked<1>();
  
  // set data directory path
  std::string data_path = "";
  std::string omp_fname = "";
  setdatapath(data_path);
  //setPdataOMP(omp_fname, 0);
  
  // prepare RNG and event
  UniformRNG rng(1);
  cgmfEvent* event = NULL;

  for ( int i = 0; i < nevents; i++ ) {
    rng.set_seed( i * seed);
    set_rng(rng);
    event = new cgmfEvent(ZAIDt, incidentEnergy, 0.0, 1.0e-8, -1);
    auto nuLF = event->getLightFragmentNeutronNu();
    auto nuHF = event->getHeavyFragmentNeutronNu();
    auto nuPre = event->getPreFissionNeutronNu();
    n[i] = nuLF + nuHF + nuPre;
    delete event;
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
