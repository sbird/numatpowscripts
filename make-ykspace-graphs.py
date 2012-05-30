import matplotlib
matplotlib.use("PDF")
import re
import sys
import numpy as np
import os.path as path
from save_figure import save_figure
import math
import matplotlib.pyplot as plt
import neutrinoplot as neu
import plot_mat_pow as pmat
import plot_flux_pow as pflux
datadir="/home/spb/data/new-kspace/yacine-kspace"
partdir="/home/spb/data/new-kspace/part"
outdir="/home/spb/papers/fs-neutrino/plots"
neut=neu.neutrino_power(data=datadir, matpow=datadir+"/b512p512nu0.3z49/")

plt.figure()

def plot_nu_powers(num_str, zz, directory):
    #High-redshift
    pmat.plot_nu_folded_power(path.join(datadir,directory+"/powerspec_nu_"+num_str+".txt"))
    pmat.plot_nu_power(path.join(datadir,directory+"/ics_matterpow_"+zz+".dat"),ls="--")
    pmat.plot_folded_power(path.join(partdir,directory+"/powerspec_type2_"+num_str+".txt"))
    #Plot delta = 1.
    xx=np.arange(1e-3,100)
    plt.loglog(xx,np.ones(np.size(xx))/xx**3,color="grey",ls="..")
    plt.xlim(1e-2,10)

plt.clf()
plot_nu_powers('008','0','b256p1024nu0.3z49')
pmat.plot_folded_power(path.join(partdir,'b256p1024nu0.3z49'+"/powerspec_type1_008.txt"))
pmat.plot_power(path.join(datadir,'b256p1024nu0.3z49'+"/ics_matterpow_0.dat"))
plt.loglog(0.2*np.ones(np.size(np.arange(1e-10,1e4))),np.arange(1e-10,1e4),color="grey",ls="--")
plt.ylim(1e-9,1000)
plt.xlabel("k /(h Mpc-1)")
plt.ylabel("P(k) (Mpc/h)$^3$")
plt.tight_layout()
save_figure(path.join(outdir,"P_nu_M03_z0"))

plt.clf()
plot_nu_powers('005','1','b256p1024nu0.3z49')
pmat.plot_folded_power(path.join(partdir,'b60p512nu0.3z49'+"/powerspec_type2_012.txt"))
plt.loglog(0.4*np.ones(np.size(np.arange(1e-12,1e4))),np.arange(1e-12,1e4),color="grey",ls="--")
plt.ylim(1e-9,1000)
plt.xlabel("k /(h Mpc-1)")
plt.ylabel("P(k) (Mpc/h)$^3$")
plt.tight_layout()
save_figure(path.join(outdir,"P_nu_M03_z1"))

plt.clf()
plot_nu_powers('003','3','b256p1024nu0.3z49')
pmat.plot_folded_power(path.join(partdir,'b60p512nu0.3z49'+"/powerspec_type2_006.txt"))
plt.loglog(np.ones(np.size(np.arange(1e-12,1e4))),np.arange(1e-12,1e4),color="grey",ls="--")
plt.ylim(1e-11,1000)
plt.xlabel("k /(h Mpc-1)")
plt.ylabel("P(k) (Mpc/h)$^3$")
plt.tight_layout()
save_figure(path.join(outdir,"P_nu_M03_z3"))

#Total matter power comparison plots
plt.clf()
#High-redshift
#Change the 512 part back to a 1024 when done.
neut.plot_directory(["/home/spb/data/new-kspace/yacine-kspace/b256p1024nu0.3z49/","/home/spb/data/new-kspace/part/b256p512nu0.3z49/"],redshifts=3,lssin=["-","-"],coloursin=["blue","red"])
plt.xlim(1e-2,10)
plt.ylim(0.7,1)
save_figure(path.join(outdir,"P_tot_M03_z3"))

#Low-redshift
plt.clf()
neut.plot_directory(["/home/spb/data/new-kspace/yacine-kspace/b256p1024nu0.3z49/","/home/spb/data/new-kspace/part/b256p512nu0.3z49/"],redshifts=0,lssin=["-","-"],coloursin=["blue","red"])
plt.xlim(1e-2,10)
plt.ylim(0.7,1)
save_figure(path.join(outdir,"P_tot_M03_z0"))


#1D Flux power
#At z=3
plt.clf()
pflux.plot_rel_flux_power(path.join(partdir,"b60p512nu0.3z49/spectra/snap_006_flux_power.txt"),path.join(datadir,"b60p512nu0z49/spectra/snap_006_flux_power.txt"),60,3,0.3,0.72,ls="--")
pflux.plot_rel_flux_power(path.join(datadir,"b60p512nu0.3z49/spectra/snap_006_flux_power.txt"),path.join(datadir,"b60p512nu0z49/spectra/snap_006_flux_power.txt"),60,3,0.3,0.72)
plt.xlim(1e-3,3e-2)
plt.ylim(1,1.05)
save_figure(path.join(outdir,"P_flux_M03_z3"))

#At z=4
plt.clf()
pflux.plot_rel_flux_power(path.join(partdir,"b60p512nu0.3z49/spectra/snap_001_flux_power.txt"),path.join(datadir,"b60p512nu0z49/spectra/snap_001_flux_power.txt"),60,4,0.3,0.72,ls="--")
pflux.plot_rel_flux_power(path.join(datadir,"b60p512nu0.3z49/spectra/snap_001_flux_power.txt"),path.join(datadir,"b60p512nu0z49/spectra/snap_001_flux_power.txt"),60,4,0.3,0.72)
plt.xlim(1e-3,3e-2)
plt.ylim(1,1.2)
save_figure(path.join(outdir,"P_flux_M03_z4"))
