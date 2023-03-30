CGMF, Cascading Gamma-ray Multiplicity and Fission
==================================================

This is a modified version of [lanl/CGMF](https://github.com/lanl/CGMF) for use in uncertainty quantification. The new capabilities include the added options listed below. Includes the following dependencies which are automatically downloaded and integrated by the build system:
- [nlohmann/json](https://github.com/nlohmann/json) for parsing json formatted parameter files
- [beykyle/omplib](https://github.com/beykyle/omplib) a library of optical model potentials and solvers

The goal of this modification is to provide an easy to use platform for running CGMF with modified internal parameters (e.g. optical model parameter, and more), and analyzing the results. 

Quickstart
-------

To build the `pyCGMF` module:

```
py setup.py build -j{nproc}
pip install -e .
```

This will:
- build `CGMF` as a library
- build the `pyCGMF` module
- build and install `CGMFtk`, if it's not already installed
- install the `pyCGMF` module

Now, we can run `CGMF`, and populate a `CGMFtk.histories.Histories` object without ever writing to disk or copying data:

```
from pyCGMF import CGMF_Input, run
from CGMFtk.histories import Histories

cgmf_input = CGMF_Input(
    nevents = 1000,
    zaid    = 98252,
    einc    = 0.0
)

histories = run(cgmf_input)
print(histories.nubar())

# we can save the histories in binary to analyze later
histories.save("histories.npy")

#...

histories.load("histories.npy")

# we can plot our results
ebins, pfns = histories.pfns()

from matplotlib import pyplot as plt

plt.step(ebins, pfns)
plt.xscale("log")
plt.xlabel(r"$E_{lab}$ [MeV]")
plt.ylabel(r"PFNS [Mev$^{-1}$]")
```

we can even run in parallel using `mpi4py`:

First, create a script called `run_cgmf.py`
```
from pyCGMF import CGMF_Input, run
from CGMFtk.histories import Histories

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# 1000 events per MPI rank
cgmf_input = CGMF_Input(
    nevents  = 1000,
    zaid     = 98252,
    einc     = 0.0,
    MPI_rank = rank
)
hist_arr = run(cgmf_input).as_numpy()
all_hist_arr = None

# gather histories from all MPI ranks to single data structure
if rank == 0:
    all_hist_arr = np.empty( (size,) + hist_arr.shape , dtype='float')
comm.Gather(hist_arr, all_hist_arr, root=0)

# create new Histories object with history info from all MPI ranks,
# write the output to disk as binary np array
all_histories = Histories(  from_arr=np.concatenate( all_histories, axis=0  ) )
all_histories.save("histories.npy")
```

This can be executed with `mpirun`:

```
mpirun -n 4 python run_cgmf.py
```

More extensive examples are contained in notebooks in `pyCGMF/examples/*.ipynb`. 

The original documentation to build CGMF as an executable is preserved below, with the addition of the extra command line arguments.

---
Authors
-------

Patrick Talou, XCP-5, Los Alamos National Laboratory, talou@lanl.gov

Ionel Stetcu, T-2, Los Alamos National Laboratory, stetcu@lanl.gov

Patrick Jaffke, Institute for Defense Analyses, pjaffke@ida.org

Michael E. Rising, XCP-3, Los Alamos National Laboratory, mrising@lanl.gov

Amy E. Lovell, T-2, Los Alamos National Laboratory, lovell@lanl.gov

Toshihiko Kawano, T-2, Los Alamos National Laboratory, kawano@lanl.gov

---
Abstract
--------


CGMF is a code that simulates the emission of prompt fission neutrons and gamma rays from excited fission fragments right after scission. It implements a Monte Carlo version of the Hauser-Feshbach statistical theory of nuclear reactions to follow the decay of the fission fragments on an event-by-event basis. Probabilities for emitting neutrons and gamma rays are computed at each stage of the decay. Each fission event history records characteristics of the parent fragment (mass, charge, kinetic energy, momentum vector, excitation energy, spin, parity) and the number (multiplicity) and characteristics (energy, direction) of the prompt neutrons and gamma rays emitted in this event.

---
Citing
------

The main publication (and documentation) to cite for CGMF is:

Patrick Talou, Ionel Stetcu, Patrick Jaffke, Michael E. Rising, Amy E. Lovell, and Toshihiko Kawano, “Fission Fragment Decay Simulations with the CGMF Code,” Comp. Phys. Comm., 269 (2021), Article 108087. DOI: [10.1016/j.cpc.2021.108087](https://doi.org/10.1016/j.cpc.2021.108087)

---
Documentation
-------------

- Patrick Talou, Ionel Stetcu, Patrick Jaffke, Michael E. Rising, Amy E. Lovell, and Toshihiko Kawano, “Fission Fragment Decay Simulations with the CGMF Code,” Comp. Phys. Comm., 269 (2021), Article 108087. Los Alamos Technical Report LA-UR-20-21264 (2020).

---
Version
-------

### Current Version 1.1.0

- Open source, BSD-3
- Copyright: Triad National Security, LLC. All rights reserved.
- Programming language: C++ (and Python for post-processing)
- Fission reactions handled: spontaneous fission of Pu-238,240,242,244 and Cf-252,254; neutron-induced fission reactions from thermal up to 20 MeV for n+U-233,234,235,238, n+Np-237, and n+Pu-239,241.


---
Building, Testing and Installing Instructions
---------------------------------------------

0) Optional, set CGMFDATA environment variable to point to data/ directory
1) Create a build directory
2) Change to build directory
3) Configuration: type `cmake ..` and then `make` for default build
    * The `..` needs to point to the top-level CGMF directory
    * Options:
        * CMAKE_BUILD_TYPE [Debug, RelWithDebInfo, Release]
        * CMAKE_INSTALL_PREFIX [path/to/install/directory]
        * cgmf.shared_library [ON/OFF (default)]
        * cgmf.x.MPI [ON/OFF (default)]
        * cgmf.tests [(default) ON/OFF]
5) Building: type `make` for default build
    * This creates the static library `libcgmf.a` in the build/libcgmf directory
    * This creates the executable `cgmf.x` in the build/utils/cgmf directory
6) Testing: type `make test` or `ctest` to run the default tests
7) Installing: type `make install` to install CGMF
    * This creates the following directory structure in the CMAKE_INSTALL_PREFIX directory:
        * bin/ [contains cgmf executable]
        * include/cgmf-1.1.0 [contains cgmf header files]
        * lib/cgmf-1.1.0 [contains libcgmf library]
        * share/data/cgmf-1.1.0 [contains cgmf data files]

Example configuration in release mode with build, test and install:
* `cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -Dcgmf.x.MPI=ON ..`
* `make && ctest && make install`
* If the build and tests pass:
    * This creates the static library `libcgmf.a` in the build directory
    * This creates the executable `cgmf.mpi.x` in the build/utils/cgmf directory
    * This installs bin, lib, cgmf/data, and cgmf/include directories into install directory

---
Execution Instructions and Options
----------------------------------

`./cgmf.x [options]`

Options:

	-i $ZAIDt   		[required]	1000*Z+A of target nucleus, or fissioning nucleus if spontaneous fission
	-e $Einc    		[required]	incident neutron energy in MeV (0.0 for spontaneous fission)
	-n $nevents 		[optional]	number of Monte Carlo fission events to run or to be read. If $nevents is negative, produces initial fission fragments yields Y(A,Z,KE,U,J,p)
	-s $startingEvent	[optional]	skip ahead to particular Monte Carlo event (1 is default)
	-f $filename		[optional]	fission histories or yields result file (default: "histories.cgmf" or "yields.cgmf")
	-t $timeCoinc		[optional]	time coincidence window for long-lived isomer gamma-ray emission cutoff (in sec)
	-d $datapath		[optional]	overrides the environment variable CGMFDATA and default datapath
	-o $omppath		[optional]	use custom global OM params, formatted in json at given path
	-r $seed                [optional]      use custom RNG seed (Note: RANDOM_SEED_BY_ARG must be set to true in config.h)
	-g $ZAIDsf              [optional]      Only de-excites fragments correspinding to ZAIDsf, with excitation energy, spin, parity, and TKE
	                                        all sampled as if that fragment was produced by the fissioning nucleus


The `-o` option currently accepts the following global optical models: [Koning-Delaroche](https://www.sciencedirect.com/science/article/pii/S0375947402013210), [WLH](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.127.182502), and [Chapel-Hill 89](https://www.sciencedirect.com/science/article/pii/037015739190039O). Example parameter files are contained in `data/optical` 

---
Results
-------

A concise summary of average results is displayed on the standard output, and event-by-event results or initial fragment distributions are saved in "histories.cgmf" or "yields.cgmf" (default names).


