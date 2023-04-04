# instructions to run:
#   mpirun -n {nproc} python -m run_cgmf_mpi
# or:
#   mpirun -n {nproc} --use-hw-threads python -m run_cgmf_mpi

from pyCGMF import CGMF_Input, run
from CGMFtk.histories import Histories

import numpy as np
from mpi4py import MPI

def run_cgmf_mpi(inp : CGMF_Input):

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    inp.MPI_rank = rank

    # run worker
    hists = run(cgmf_input)

    # write results from just this rank:
    hist_arr.save("histories_rank_{}.npy".format(rank))

    # gather histories from all MPI ranks
    result = comm.gather(hist_arr.histories, root=0)

    # concatenate them and print the result
    if rank == 0:
        # create new Histories object with history info from all MPI ranks,
        # write the output to disk as binary np array
        all_histories = Histories(  from_arr=np.concatenate( result, axis=0  ) )
        all_histories.save("all_histories.npy")
        return all_histories

    return None

def test():
    for i in range(size):

        # just check columns 0-4: A, Z, J, P, U
        test = np.array(
            Histories.load("histories_rank_{}.npy".format(i)).histories[:,0:5]
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

def main():
    # 100 events per MPI rank
    cgmf_input = CGMF_Input(
        nevents  = 100,
        zaid     = 98252,
        einc     = 0.0
    )

    # run on all workers
    # result will only be meaningful on rank 0
    result = run_cgmf_mpi(inp)

    if rank == 0:
        print("nubar: {}".format(result.nubartot()))
        test()

if __name__ == "__main__":
    main()
