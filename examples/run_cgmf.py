from pyCGMF import Input, run
from CGMFtk.histories import Histories
import numpy as np

cgmf_input = Input(
    nevents = 100,
    zaid    = 98252,
    einc    = 0.0
)

histories = run(cgmf_input)
print("nubar: {}".format( histories.nubartot() ) )

# we can save the histories to disk (in compact binary fmt using numpy.save/load)
histories.save("histories.npy")

#...
# later, we can load them back into memory
hist2 = Histories.load("histories.npy")

# we can plot our results
ebins, pfns = hist2.pfns()

from matplotlib import pyplot as plt

plt.step(ebins, pfns)
plt.xscale("log")
plt.xlabel(r"$E_{lab}$ [MeV]")
plt.ylabel(r"PFNS [Mev$^{-1}$]")
plt.show()
