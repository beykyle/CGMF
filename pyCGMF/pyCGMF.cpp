#ifndef PYBINDINGS
#define PYBINDINGS

#include <string>

#include <pybind11/pybind11.h>
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

template<typename T>
using arr_t  = typename py::array_t<T>;
template<typename T>
using shape_t = typename py::array_t<T>::ShapeContainer;

struct CGMF_Input {
  int nevents;
  int ZAIDt;
  double einc;
  
  double time_coinc_wndw   = 1.0e-8;
  int MPI_rank             = 0;
  int seed                 = 13;
  std::string omp_fpath    = "";
};

struct EventData {
  arr_t<double>   fragments;  
  arr_t<py::list> neutrons;   
  arr_t<py::list> gammas;
  arr_t<py::list> pre_fission_neutrons;

  EventData(const CGMF_Input& input) {
    // 2 fragments / event * n events
    const auto N = input.nevents * 2;
    
    // fragments has 16 fields in CGMFtk.histories.Histories (A,Z,J, ...)
    fragments.resize( shape_t<double>{ N, 16 } );
    
    // fields for emitted particles - size depends on event multiplicity;
    // hence the list.
    // neutrons have 4 fields: e_cm, e_lab, dircos_cm, dircos_lab
    const auto nshape   = shape_t<py::list>{ N, 4 };
    // gammas have e_cm, e_lab, ages
    const auto gshape   = shape_t<py::list>{ N, 3 };
    // pre fission neutrons have 2 fields: e_lab, dircos_lab
    const auto pfnshape = shape_t<py::list>{ N, 2 };
    
    neutrons.resize(nshape);
    gammas.resize(gshape);            
    pre_fission_neutrons.resize(pfnshape);
  }

  void process(cgmfEvent* event, int i) {
    
    auto frags = fragments.mutable_unchecked<2>();
    auto ns    = neutrons.mutable_unchecked<2>();
    auto gs    = neutrons.mutable_unchecked<2>();
    auto pfn   = pre_fission_neutrons.mutable_unchecked<2>();

    auto idx = 2 * i;
    
    const auto  nul   = event->getLightFragmentNeutronNu();
    const auto  nugl  = event->getLightFragmentPhotonNu();
    const auto  nuh   = event->getHeavyFragmentNeutronNu();
    const auto  nugh  = event->getHeavyFragmentPhotonNu();
    const auto  nup   = event->getPreFissionNeutronNu();

    // light fragment
    frags(idx,0) = event->getLightFragmentMass();
    frags(idx,1) = event->getLightFragmentCharge();
    frags(idx,2) = event->getLightFragmentExcitationEnergy();
    frags(idx,3) = event->getLightFragmentSpin();
    frags(idx,4) = event->getLightFragmentParity();
    frags(idx,5) = event->getLightFragmentKineticEnergy();
    frags(idx,6) = nul;
    frags(idx,7) = nugl;

    frags(idx,8)  = event->getLightFragmentPreMomentumX();
    frags(idx,9)  = event->getLightFragmentPreMomentumY();
    frags(idx,10) = event->getLightFragmentPreMomentumZ();
    frags(idx,11) = event->getLightFragmentPostMomentumX();
    frags(idx,12) = event->getLightFragmentPostMomentumY();
    frags(idx,13) = event->getLightFragmentPostMomentumZ();
    frags(idx,14) = 0; 
    frags(idx,15) = event->getLightFragmentKineticEnergyPost();
    
    // neutrons LF
    for (int n = 0; n < nul; ++ n) {
      arr_t<double> dircos_cm (
          std::array<double,3>{
            event->getCmNeutronDircosu(n),
            event->getCmNeutronDircosv(n),
            event->getCmNeutronDircosw(n)
          }
        );
      arr_t<double> dircos_lab (
          std::array<double,3>{
            event->getNeutronDircosu(n),
            event->getNeutronDircosv(n),
            event->getNeutronDircosw(n)
          }
        );

      ns(idx,0).append(event->getCmNeutronEnergy(n));
      ns(idx,1).append(event->getNeutronEnergy(n));
      ns(idx,2).append(dircos_cm);
      ns(idx,3).append(dircos_lab);
    }
    
    // gammas LF
    for (int n = 0; n < nugl; ++ n) {
      // TODO ignore doppler shift?
      ns(idx,0).append(event->getPhotonEnergy(n));
      ns(idx,1).append(event->getPhotonEnergy(n));
      ns(idx,3).append(event->getPhotonAge(n));
    }

    idx++;
    
    // heavy fragment
    frags(idx,0) = event->getHeavyFragmentMass();
    frags(idx,1) = event->getHeavyFragmentCharge();
    frags(idx,2) = event->getHeavyFragmentExcitationEnergy();
    frags(idx,3) = event->getHeavyFragmentSpin();
    frags(idx,4) = event->getHeavyFragmentParity();
    frags(idx,5) = event->getHeavyFragmentKineticEnergy();
    frags(idx,6) = nuh;
    frags(idx,7) = nugh;

    frags(idx,8)  = event->getHeavyFragmentPreMomentumX();
    frags(idx,9)  = event->getHeavyFragmentPreMomentumY();
    frags(idx,10) = event->getHeavyFragmentPreMomentumZ();
    frags(idx,11) = event->getHeavyFragmentPostMomentumX();
    frags(idx,12) = event->getHeavyFragmentPostMomentumY();
    frags(idx,13) = event->getHeavyFragmentPostMomentumZ();
    frags(idx,14) = nup;
    frags(idx,15) = event->getHeavyFragmentKineticEnergyPost();

    // neutrons HF
    for (int n = nul; n < nul + nuh; ++ n) {
      arr_t<double> dircos_cm (
          std::array<double,3>{
            event->getCmNeutronDircosu(n),
            event->getCmNeutronDircosv(n),
            event->getCmNeutronDircosw(n)
          }
        );
      arr_t<double> dircos_lab (
          std::array<double,3>{
            event->getNeutronDircosu(n),
            event->getNeutronDircosv(n),
            event->getNeutronDircosw(n)
          }
        );

      ns(idx,0).append(event->getCmNeutronEnergy(n));
      ns(idx,1).append(event->getNeutronEnergy(n));
      ns(idx,2).append(dircos_cm);
      ns(idx,3).append(dircos_lab);
    }
    
    // gammas HF
    for (int n = nugl; n < nugl + nugh; ++ n) {
      // TODO ignore doppler shift?
      ns(idx,0).append(event->getPhotonEnergy(n));
      ns(idx,1).append(event->getPhotonEnergy(n));
      ns(idx,3).append(event->getPhotonAge(n));
    }
    
    // pre-fission neutrons
    for (int n = 0; n < nup; ++ n) {
      arr_t<double> dircos_lab (
          std::array<double,3>{
            event->getPreFissionNeutronDircosu(n),
            event->getPreFissionNeutronDircosv(n),
            event->getPreFissionNeutronDircosw(n)
          }
        );

      pfn(idx,0).append(event->getPreFissionNeutronEnergy(n));
      pfn(idx,1).append(dircos_lab);
    }
  }
};

EventData run(const CGMF_Input& input) {

  // pre-allocate data arrays based on number of events
  auto event_data = EventData(input);
  
  // set data directory path
  std::string data_path = "";
  setdatapath(data_path);
  
  // check if opical model param file path has been passes
  if (not input.omp_fpath.empty())
    setPdataOMP(input.omp_fpath, input.MPI_rank);
  
  // initialize RNG and event
  UniformRNG rng(1);
  cgmfEvent* event = NULL;

  for ( int i = 0; i < input.nevents; i++ ) {
    rng.set_seed( input.seed + i + i * input.seed);
    set_rng(rng);
    event = new cgmfEvent(
        input.ZAIDt, input.einc, 0.0, input.time_coinc_wndw, -1);
    
    // fill data for this event
    event_data.process(event, i);
    
    delete event;
  }
  return event_data;
}


PYBIND11_MODULE(pyCGMF, m) {
  m.doc() = "pyCGMF: python bindings for running CGMF on MPI";
  
  m.def(
      "run", 
      &run, 
      "runs nevents CGMF histories of 252-Cf (sf), returns data relating to each event" 
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
