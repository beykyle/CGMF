#ifndef PYBINDINGS
#define PYBINDINGS

#include <string>

#include <pybind11/cast.h>
#include <pybind11/pytypes.h>
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

using namespace pybind11::literals;
namespace py = pybind11;

py::module_ np = py::module_::import("numpy"); 
py::module_ CGMFtk = py::module_::import("CGMFtk");

template<typename T>
using arr_t  = typename py::array_t<T>;
template<typename T>
using shape_t = typename py::array_t<T>::ShapeContainer;

struct CGMF_Input {
  int nevents;
  int zaid;
  double einc;
  
  double time_coinc_wndw   = 1.0e-8;
  int MPI_rank             = 0;
  int seed                 = 13;
  std::string omp_fpath    = "";
};

enum class Fragment : int {
  A, Z, U, J, P, KEPre, Nu, Nug, PxPre, PyPre, PzPre, PxPost, PyPost, PzPost, NuPFN, KEPost
};

struct EventData {
  //CGMFtk data that us fixed in size we can store as np array
  arr_t<double> fragments;

  // CGMFtk data that is not fixed in size (e.g. neutron energies depends on 
  // CGMF's simulated neutron multiplicities), we store as lists of lists
  py::list neutron_Elab;
  py::list neutron_Ecm;
  py::list gamma_Elab;
  py::list gamma_Ecm;
  py::list gamma_Age;
  py::list neutron_dircoslab;
  py::list neutron_dircoscm;
  py::list pf_neutron_Elab;
  py::list pf_neutron_dircoslab;

  EventData(const CGMF_Input& input) {
    // 2 fragments / event * n events
    const auto N = input.nevents * 2;
    
    // 16 fields enumerated in Fragment
    fragments.resize( shape_t<double>{ 16, N } );
  }

  py::array get_field_data( Fragment field ) {
    const auto i = static_cast<int>(field);
    auto res = py::array{ fragments[py::make_tuple(i, py::ellipsis())] };
    return np.attr("array")(res);
  }

  void process(cgmfEvent* event, int i) {
    
    auto frags = fragments.mutable_unchecked<2>();
    
    // initialize inner lists for emitted particle data in event
    // neutron
    auto comEn = py::list();
    auto labEn = py::list();
    auto labdc = py::list();
    auto comdc = py::list();
    // gamma
    auto comEg = py::list();
    auto labEg = py::list();
    auto gAge  = py::list();
    // pre-fission neutron
    auto pflabEn = py::list();
    auto pflabdc = py::list();

    auto idx = 2 * i;
    
    const auto  nul   = event->getLightFragmentNeutronNu();
    const auto  nugl  = event->getLightFragmentPhotonNu();
    const auto  nuh   = event->getHeavyFragmentNeutronNu();
    const auto  nugh  = event->getHeavyFragmentPhotonNu();
    const auto  nup   = event->getPreFissionNeutronNu();

    // light fragment
    frags(Fragment::A,idx)     = event->getLightFragmentMass();
    frags(Fragment::Z,idx)     = event->getLightFragmentCharge();
    frags(Fragment::U,idx)     = event->getLightFragmentExcitationEnergy();
    frags(Fragment::J,idx)     = event->getLightFragmentSpin();
    frags(Fragment::P,idx)     = event->getLightFragmentParity();
    frags(Fragment::KEPre,idx) = event->getLightFragmentKineticEnergy();
    frags(Fragment::Nu,idx)    = nul;
    frags(Fragment::Nug,idx)   = nugl;

    frags(Fragment::PxPre,idx)  = event->getLightFragmentPreMomentumX();
    frags(Fragment::PyPre,idx)  = event->getLightFragmentPreMomentumY();
    frags(Fragment::PzPre,idx)  = event->getLightFragmentPreMomentumZ();
    frags(Fragment::PxPost,idx) = event->getLightFragmentPostMomentumX();
    frags(Fragment::PyPost,idx) = event->getLightFragmentPostMomentumY();
    frags(Fragment::PzPost,idx) = event->getLightFragmentPostMomentumZ();
    frags(Fragment::NuPFN,idx)  = 0;
    frags(Fragment::KEPost,idx) = event->getLightFragmentKineticEnergyPost();
    
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
      
      comEn.append( event->getCmNeutronEnergy(n) );
      labEn.append( event->getNeutronEnergy(n) );
      labdc.append( dircos_lab);
      comdc.append( dircos_cm);
    }
    
    // gammas LF
    for (int n = 0; n < nugl; ++ n) {
      // TODO do we ignore doppler shift in CGMF?
      comEg.append( event->getPhotonEnergy(n) );
      labEg.append( event->getPhotonEnergy(n) );
      gAge.append(  event->getPhotonAge(n) );
    }
    
    // append our inner emitted particle lists for the LF to the event-by-event lists
    neutron_Elab.append( labEn );
    neutron_Ecm.append( comEn );
    gamma_Elab.append( labEg );
    gamma_Ecm.append( comEg );
    gamma_Age.append( gAge );
    neutron_dircoslab.append( labdc );
    neutron_dircoscm.append( comdc );
    pf_neutron_Elab.append( pflabEn );
    pf_neutron_dircoslab.append( pflabdc );

    // clear inner lists
    comEn = py::list();
    labEn = py::list();
    labdc = py::list();
    comdc = py::list();
    comEg = py::list();
    labEg = py::list();
    gAge  = py::list();
    pflabEn = py::list();
    pflabdc = py::list();

    idx++;
    
    // heavy fragment
    frags(Fragment::A,idx)     = event->getHeavyFragmentMass();
    frags(Fragment::Z,idx)     = event->getHeavyFragmentCharge();
    frags(Fragment::U,idx)     = event->getHeavyFragmentExcitationEnergy();
    frags(Fragment::J,idx)     = event->getHeavyFragmentSpin();
    frags(Fragment::P,idx)     = event->getHeavyFragmentParity();
    frags(Fragment::KEPre,idx) = event->getHeavyFragmentKineticEnergy();
    frags(Fragment::Nu,idx)    = nuh;
    frags(Fragment::Nug,idx)   = nugh;

    frags(Fragment::PxPre,idx)  = event->getHeavyFragmentPreMomentumX();
    frags(Fragment::PyPre,idx)  = event->getHeavyFragmentPreMomentumY();
    frags(Fragment::PzPre,idx)  = event->getHeavyFragmentPreMomentumZ();
    frags(Fragment::PxPost,idx) = event->getHeavyFragmentPostMomentumX();
    frags(Fragment::PyPost,idx) = event->getHeavyFragmentPostMomentumY();
    frags(Fragment::PzPost,idx) = event->getHeavyFragmentPostMomentumZ();
    frags(Fragment::NuPFN,idx)  = nup;
    frags(Fragment::KEPost,idx) = event->getHeavyFragmentKineticEnergyPost();

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
      
      comEn.append( event->getCmNeutronEnergy(n) );
      labEn.append( event->getNeutronEnergy(n) );
      labdc.append( dircos_lab);
      comdc.append( dircos_cm);

    }
    
    // gammas HF
    for (int n = nugl; n < nugl + nugh; ++ n) {
      comEg.append( event->getPhotonEnergy(n) );
      labEg.append( event->getPhotonEnergy(n) );
      gAge.append(  event->getPhotonAge(n) );
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
      
      pflabEn.append( event->getPreFissionNeutronEnergy(n) );
      pflabdc.append( dircos_lab );
    }
    if (nup == 0) {
      // np array w/ dtype 'object' filled with all empty lists
      // misbehaves, so we have to fill it with garbage if
      // there are no pre-fission neutrons
      pflabEn.append( 0.0 );
      pflabdc.append( arr_t<double>({0.,0.,0.}));
    }
    
    // append our inner emitted particle lists for the HF to the event-by-event lists
    neutron_Elab.append( labEn );
    neutron_Ecm.append( comEn );
    gamma_Elab.append( labEg );
    gamma_Ecm.append( comEg );
    gamma_Age.append( gAge );
    neutron_dircoslab.append( labdc );
    neutron_dircoscm.append( comdc );
      pf_neutron_Elab.append( pflabEn );
      pf_neutron_dircoslab.append( pflabdc );
  }

  py::object concat() {
    // slap that numpy array in the format CGMFtk expects it
    return np.attr("dstack")(
        py::make_tuple(
          get_field_data(Fragment::A), 
          get_field_data(Fragment::Z),
          get_field_data(Fragment::U),
          get_field_data(Fragment::J),
          get_field_data(Fragment::P),
          get_field_data(Fragment::KEPre),
          get_field_data(Fragment::Nu),
          get_field_data(Fragment::Nug),
          neutron_Ecm,
          neutron_Elab,
          gamma_Ecm,
          gamma_Elab,
          gamma_Age,
          get_field_data(Fragment::PxPre),
          get_field_data(Fragment::PyPre),
          get_field_data(Fragment::PzPre),
          get_field_data(Fragment::PxPost),
          get_field_data(Fragment::PyPost),
          get_field_data(Fragment::PzPost),
          neutron_dircoscm,
          neutron_dircoslab,
          get_field_data( Fragment::NuPFN),
          pf_neutron_Elab,
          pf_neutron_dircoslab,
          get_field_data(Fragment::KEPost) 
        ))[py::make_tuple(0, py::ellipsis())];
  }

};

py::object run(const CGMF_Input& input) {

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
        input.zaid, input.einc, 0.0, input.time_coinc_wndw, -1);
    
    // fill data for this event
    event_data.process(event, i);
    
    delete event;
  }
  return CGMFtk.attr("Histories")(
      "from_arr"_a =  event_data.concat() );
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
       py::arg("zaid") = 98252,
       py::arg("einc") = 0.0,
       py::arg("time_coinc_wndw") = 1.0e-8,
       py::arg("MPI_rank") = 0,
       py::arg("seed") = 13,
       py::arg("omp_fpath") = ""
      )
  .def_readwrite("nevents", &CGMF_Input::nevents)
  .def_readwrite("zaid", &CGMF_Input::zaid)
  .def_readwrite("einc", &CGMF_Input::einc)
  .def_readwrite("time_coinc_wndw", &CGMF_Input::time_coinc_wndw)
  .def_readwrite("MPI_rank", &CGMF_Input::MPI_rank)
  .def_readwrite("seed", &CGMF_Input::seed)
  .def_readwrite("omp_fpath", &CGMF_Input::omp_fpath);
}

#endif
