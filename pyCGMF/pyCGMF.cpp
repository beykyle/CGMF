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

namespace py = pybind11;
using arr = py::array_t<double>;

struct CGMF_Input {
  int nevents;
  int ZAIDt;
  double einc;
  
  double time_coinc_wndw   = 1.0e-8;
  int MPI_rank             = 0;
  int seed                 = 13;
  std::string omp_fpath    = "";
};

arr run_histories_get_nubar(const CGMF_Input& input) {
  // prepare array
  auto nu = arr();
  nu.resize( arr::ShapeContainer{input.nevents} );
  auto n = nu.mutable_unchecked<1>();
  
  // set data directory path
  std::string data_path = "";
  std::string omp_fpath = "";
  setdatapath(data_path);
  if (not input.omp_fpath.empty())
    setPdataOMP(input.omp_fpath, input.MPI_rank);
  
  // prepare RNG and event
  UniformRNG rng(1);
  cgmfEvent* event = NULL;

  for ( int i = 0; i < input.nevents; i++ ) {
    rng.set_seed( input.seed + i + i * input.seed);
    set_rng(rng);
    event = new cgmfEvent(
        input.ZAIDt, input.einc, 0.0, input.time_coinc_wndw, -1);
    auto nuLF = event->getLightFragmentNeutronNu();
    auto nuHF = event->getHeavyFragmentNeutronNu();
    auto nuPre = event->getPreFissionNeutronNu();
    n[i] = nuLF + nuHF + nuPre;
    delete event;
  }
  return nu;
}


PYBIND11_MODULE(pyCGMF, m) {
  m.doc() = "pyCGMF: python bindings for running CGMF on MPI";
  
  m.def(
      "run_histories_get_nubar", 
      &run_histories_get_nubar, 
      "runs nevents CGMF histories of 252-Cf (sf), returns the neutron multiplicity for each event" 
    );
  
  py::class_<CGMF_Input>(m, "CGMF_Input")
  .def(
       py::init<int, int, double, double, int, int, std::string >(),
       py::arg("nevents") = 100,
       py::arg("ZAIDt") = 98252,
       py::arg("einc") = 0.0,
       py::arg("time_coinc_wndw") = 1.0e-8,
       py::arg("MPI_rank") = 0,
       py::arg("seed") = 13,
       py::arg("omp_fpath") = ""
      )
  .def_readwrite("nevents", &CGMF_Input::nevents)
  .def_readwrite("ZAID_t", &CGMF_Input::ZAIDt)
  .def_readwrite("einc", &CGMF_Input::einc)
  .def_readwrite("time_coinc_wndw", &CGMF_Input::time_coinc_wndw)
  .def_readwrite("MPI_rank", &CGMF_Input::MPI_rank)
  .def_readwrite("seed", &CGMF_Input::seed)
  .def_readwrite("omp_fpath", &CGMF_Input::omp_fpath);
}

#endif
