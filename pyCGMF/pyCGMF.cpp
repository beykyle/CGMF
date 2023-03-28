#ifndef PYBINDINGS
#define PYBINDINGS

#include <iostream>

#include <pybind11/pybind11.h>

#include "solver/scatter.hpp"


namespace py = pybind11;

PYBIND11_MODULE(pyCGMF, m) {
  m.doc() = "pyCGMF: python bindings for running CGMF on MPI";
  m.def("test", [](){ return "test"; }, "returns the string \"test\"" );
}

#endif
