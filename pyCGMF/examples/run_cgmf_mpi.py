# instructions to run:
#   mpirun -n {nproc} python -m run_cgmf_mpi
# or:
#   mpirun -n {nproc} --use-hw-threads python -m run_cgmf_mpi

from pyCGMF import CGMF_Input, run
from CGMFtk.histories import Histories

import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# 1000 events per MPI rank
cgmf_input = CGMF_Input(
    nevents  = 100,
    zaid     = 98252,
    einc     = 0.0,
    MPI_rank = rank
)

hist_arr = run(cgmf_input)

# write results from just this rank:
np.save( "histories_rank_{}.npy".format(rank), hist_arr, allow_pickle=True)

# gather histories from all MPI ranks to single data structure
result = comm.gather(hist_arr, root=0)

if rank == 0:

    # create new Histories object with history info from all MPI ranks,
    # write the output to disk as binary np array
    all_histories = Histories(  from_arr=np.concatenate( result, axis=0  ) )
    all_histories.save("all_histories.npy")

    # check to make sure histories are all where they should be
    for i in range(size):

        # just check columns 0-4: A, Z, J, P, U
        test = np.array(
            np.load("histories_rank_{}.npy".format(i), allow_pickle=True)[:,0:5]
          , dtype="float")

        # MPI gather respects the ordering of the ranks, so we know
        # exactly what indices each rank's results will be in after
        # concatenation
        ibot = i * cgmf_input.nevents*2
        itop = ibot + cgmf_input.nevents*2
        gathered =np.array(
            all_histories.histories[ibot:itop,0:5],
            dtype="float")

        # comparison
        assert( np.allclose(test, gathered))


