import matplotlib
matplotlib.use("PDF")
import re
import sys
import numpy as np
import os.path as path
from save_figure import save_figure
import math
import h5py
import matplotlib.pyplot as plt
import plot_mat_pow as pmat
datadir="/home/spb/data/hybrid-kspace/"
outdir="/home/spb/data/hybrid-kspace/plots"

plt.figure()

def plot_snap(num):
    pmat.plot_nu_folded_power(path.join(datadir,"b512p512nu0.45z99/output/powerspec_nu_"+num+".txt"),ls="--",color="black", label="Fourier-space")
    pmat.plot_folded_power(path.join(datadir,"b512p512nu0.45z99/outputnu/powerspec_type2_"+num+".txt"),ls="-",color="blue", label="Hybrid")
    pmat.plot_folded_power(path.join(datadir,"b512p512nu0.45z99/outputnu3/powerspec_type2_"+num+".txt"),ls="-",color="red", label="Hybridz3")
    pmat.plot_folded_power(path.join(datadir,"b512p512nu0.45z99part/output/powerspec_type2_"+num+".txt"),ls="--",color="grey", label="Particle")
    plt.legend()
    save_figure(path.join(outdir,"hyb_part_n0.45_"+num))
    plt.clf()

plot_snap('005')
plot_snap('008')
plot_snap('007')
plot_snap('006')

def vel_hist(sdir, num, bins=80, label="None"):
    vels = []
    for x in xrange(500):
        try:
            f = h5py.File(path.join(sdir, "snapdir_"+num+"/snap_"+num+"."+str(x)+".hdf5"),'r')
            vels.append(np.sqrt(np.sum(np.array(f["PartType2"]["Velocities"])**2,axis=1)))
            f.close()
        except IOError:
            break
    vel = np.concatenate(vels)
    (nn, binsout) = np.histogram(np.log10(vel), bins=np.log10(bins))
    plt.semilogx(bins[1:], nn, '--', label=label)
    return (nn, binsout)

def plot_v_hist(num):
    plt.figure(1)
    bins = np.logspace(1, 4.6,140)
    (nn1, _) = vel_hist(path.join(datadir,"b512p512nu0.45z99/outputnu"), num,bins=bins, label="Hybrid")
    (nn2, _) = vel_hist(path.join(datadir,"b512p512nu0.45z99part/output"), num,bins=bins, label="Part")
    plt.legend(loc="upper left")
    plt.xlim(100, 10**4)
    save_figure(path.join(outdir, "velhist_n0.45_z"+num))
    plt.clf()
    plt.figure(9)
    plt.semilogx(bins[1:],1.*nn1/(nn2+0.001), '-', label=num)

plt.figure(1)
bins = np.logspace(1, 4.6,140)
# (nn1, _) = vel_hist(path.join(datadir,"b512p512nu0.45z99/outputnu3"), '005',bins=bins, label="Hybrid z3")
(nn2, _) = vel_hist(path.join(datadir,"b512p512nu0.45z99/outputnunl"), '005',bins=bins, label="Hybrid nonlinear")
plot_v_hist('005')
plot_v_hist('006')
plot_v_hist('007')
plot_v_hist('008')
plt.figure(9)
plt.legend()
plt.ylim(0.5,1.5)
save_figure(path.join(outdir, "veldiff_n0.45"))
plt.clf()
# vel_hist(path.join(datadir,"b512p512nu0.45z99/outputnu"), '000',"Hybrid")
# vel_hist(path.join(datadir,"b512p512nu0.45z99part/output"), '000',"Part")
# plt.legend(loc="upper left")
# plt.xlim(1e4, 1e8)
# save_figure(path.join(outdir, "velhist_n0.45_z000"))
# plt.clf()
pmat.plot_folded_power(path.join(datadir,'b512p512nu0.45z99/output/powerspec_008.txt'),ls="-", color="black")
pmat.plot_power(path.join(datadir,"b512p512nu0.45z99/ics_matterpow_0.dat"))
save_figure(path.join(outdir, "kspace_linear_n0.45_z0"))
plt.clf()
pmat.plot_nu_folded_power(path.join(datadir,"b512p512nu0.45z99/output/powerspec_nu_005.txt"),ls="--",color="black", label="Fourier-space")
pmat.plot_genpk_single_power(path.join(datadir,"b512p512nu0.45z99/outputnu/snapdir_005/PK-nu-snap_005"),512,ls="-",color="blue", label="Hybrid")
# pmat.plot_folded_power(path.join(datadir,"b512p512nu0.45z99/outputnu/powerspec_type2_005.txt"),ls="--",color="blue", label="Hybrid")
pmat.plot_folded_power(path.join(datadir,"b512p512nu0.45z99part/output/powerspec_type2_005.txt"),ls="--",color="grey", label="Particle")
plt.legend()
save_figure(path.join(outdir,"hyb_part_n0.45_z1"))
plt.clf()
pmat.plot_folded_power(path.join(datadir,'b512p512nu0.45z99/output/powerspec_005.txt'),ls="-", color="black")
pmat.plot_power(path.join(datadir,"b512p512nu0.45z99/ics_matterpow_1.dat"))
save_figure(path.join(outdir, "kspace_linear_n0.45_z1"))
plt.clf()
