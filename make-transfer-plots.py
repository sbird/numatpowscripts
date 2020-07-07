import matplotlib
matplotlib.use("PDF")
import re
import sys
import numpy as np
import os.path as path
from save_figure import save_figure
import math
import matplotlib.pyplot as plt
# import neutrinoplot as neu
import plot_mat_pow as pmat
# import plot_flux_pow as pflux
outdir="/home/spb/data/hybrid-kspace/plots"
nudatadir = "/home/spb/data/hybrid-kspace/test2/b512p512nu0.45z99hyb/camb_linear"
dmdatadir = "/home/spb/data/hybrid-kspace/test2/b512p512nu0z99/camb_linear"

(knu, pknu_z2) = pmat.get_power(path.join(nudatadir,"ics_matterpow_2.dat"))
(kdm, pkdm_z2) = pmat.get_power(path.join(dmdatadir,"ics_matterpow_2.dat"))
pkdm_z2 = pmat.rebin(pkdm_z2, kdm, knu)
plt.semilogx(knu, pknu_z2/pkdm_z2, ls="-", color="blue", label="z=2")
(knu, pknu_z0) = pmat.get_power(path.join(nudatadir,"ics_matterpow_0.dat"))
(kdm, pkdm_z0) = pmat.get_power(path.join(dmdatadir,"ics_matterpow_0.dat"))
pkdm_z0 = pmat.rebin(pkdm_z0, kdm, knu)
plt.semilogx(knu, pknu_z0/pkdm_z0, ls="--", color="black", label="z=0")
plt.xlim(1e-3, 10)
plt.legend()
save_figure(path.join(outdir, "transfer_redshift"))
# plot_genpk_rel_power(matpow1,matpow2, box,o_nu = 0, colour="blue"):
