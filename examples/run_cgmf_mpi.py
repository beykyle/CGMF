# instructions to run:
#   mpirun -n {nproc} python -m run_cgmf_mpi.py
# or:
#   mpirun -n {nproc} --use-hw-threads python -m run_cgmf_mpi.py

from pyCGMF import CGMF_Input, run
from CGMFtk.histories import Histories

import sys
import numpy as np
from mpi4py import MPI


def run_cgmf_mpi(inp: CGMF_Input, comm):
    print("Running {} histories on rank {}".format(inp.nevents, inp.MPI_rank))
    sys.stdout.flush()
    # run worker
    hists = run(inp)

    # write results from just this rank:
    print("Printing results from rank {}".format(inp.MPI_rank))
    sys.stdout.flush()
    hists.save("histories_rank_{}.npy".format(inp.MPI_rank))

    # gather histories from all MPI ranks
    result = comm.gather(hists.histories, root=0)

    # concatenate them and print the result
    if inp.MPI_rank == 0:
        # create new Histories object with history info from all MPI ranks,
        # write the output to disk as binary np array
        all_histories = Histories(from_arr=np.concatenate(result, axis=0))
        all_histories.save("all_histories.npy")
        return all_histories

    return None


def test(inp: CGMF_Input, size: int, all_histories: Histories):
    for i in range(size):
        # just check columns 0-4: A, Z, J, P, U
        test = np.array(
            Histories.load("histories_rank_{}.npy".format(i)).histories[:, 0:5],
            dtype="float",
        )

        # MPI gather respects the ordering of the ranks, so we know
        # exactly what indices each rank's results will be in after
        # concatenation
        ibot = i * inp.nevents * 2
        itop = ibot + inp.nevents * 2
        gathered = np.array(all_histories.histories[ibot:itop, 0:5], dtype="float")

        # comparison
        assert np.allclose(test, gathered)


def main():
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    # 100 events per MPI rank
    inp = CGMF_Input(nevents=100, zaid=98252, einc=0.0)

    inp.MPI_rank = rank

    # run on all workers
    result = run_cgmf_mpi(inp, comm)

    if rank == 0:
        assert result is not None  # should retun Histories instance on rank 0
        test(inp, size, result)
        print("Test passed!")
        print("nubar: {}".format(result.nubartot()))


if __name__ == "__main__":
    main()
