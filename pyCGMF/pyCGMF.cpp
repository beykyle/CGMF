#ifndef PYBINDINGS
#define PYBINDINGS

#include <string>

#include <pybind11/pybind11.h>

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

//namespace py = pybind11;

PYBIND11_MODULE(pyCGMF, m) {
  m.doc() = "pyCGMF: python bindings for running CGMF on MPI";
  m.def("test", &test_func, "returns the string \"pyCGMF test\"" );
}

#endif
