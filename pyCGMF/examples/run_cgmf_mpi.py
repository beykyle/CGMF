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
