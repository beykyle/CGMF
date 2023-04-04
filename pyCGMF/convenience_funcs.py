from pyCGMF import CGMF_Input, run
from CGMFtk.histories import Histories

import numpy as np
from mpi4py import MPI


def run_cgmf(inp : CGMF_Input):
    """
    Run CGMF and get the resulting CGMFtk.histories.Histories object
    """
    return Histories( from_arr=run(inp) )

def run_cgmf_mpi(inp : CGMF_Input):
    """
    Run CGMF in parallel with MPI and get the resulting CGMFtk.histories.Histories object
    """
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    inp.MPI_rank = rank

    # run on worker
    hists = run(inp)

    # concatenate and return results from main
    result = comm.gather(hists, root=0)
    if rank == 0:
        return Histories(  from_arr=np.concatenate( result, axis=0  ) )
