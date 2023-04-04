from pyCGMF import CGMF_Input, run
from CGMFtk.histories import Histories
import numpy as np

cgmf_input = CGMF_Input(
    nevents = 100,
    zaid    = 98252,
    einc    = 0.0
)

histories = Histories( from_arr=run(cgmf_input) )
print("nubar: {}".format( histories.nubar() ) )

# we can save the histories in binary to analyze later
histories.save("histories.npy")

#...

histories = Histories.load("histories.npy")

# we can plot our results
ebins, pfns = histories.pfns()

from matplotlib import pyplot as plt

plt.step(ebins, pfns)
plt.xscale("log")
plt.xlabel(r"$E_{lab}$ [MeV]")
plt.ylabel(r"PFNS [Mev$^{-1}$]")
plt.show()
