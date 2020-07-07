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
datadir="/home/spb/data/hybrid-kspace/test2"
ykdatadir="/home/spb/data2/new-kspace"
outdir="/home/spb/data/hybrid-kspace/plots"

plt.figure()

def plot_nu_powers(datadir, num_str, zz, directory, dir2):
    """Plot the neutrino power spectrum from both semi-linear and particle codes."""
    pmat.plot_nu_folded_power(path.join(datadir,dir2+"/output/powerspec_nu_"+num_str+".txt"),ls="-",color="black")
    (k,delta)=pmat.get_folded_power(path.join(datadir,directory+"/output/powerspec_type2_"+num_str+".txt"))
    box=512
    part=512
    shot=(box/part)**3/(2*math.pi**2)*np.ones(np.size(delta))
#     plt.loglog(k,shot,color="grey",ls=":")
    plt.loglog(k,delta,ls="--",color="pink")
    plt.loglog(k,delta-shot,ls="--",color="red")
    try:
        pmat.plot_nu_power(path.join(datadir,dir2+"/camb_linear/ics_matterpow_"+zz+".dat"),ls=":")
    except Exception:
        pass
#     plt.loglog(k,delta-(256./shot)**3/(2*math.pi**2)*np.ones(np.size(delta)),ls="--")
    plt.xlim(1e-2,10)

def plot_hyb_nu_power(num_str, datadir, directory, part_prop=0.1):
    (k_part,pk_part)=pmat.get_folded_power(path.join(datadir,directory+"/output/powerspec_type2_"+num_str+".txt"))
    (k_sl, pk_sl) = pmat.get_nu_folded_power(path.join(datadir,directory+"/output/powerspec_nu_"+num_str+".txt"))
    pk_part_r = pmat.rebin(pk_part, k_part, k_sl)
    shot=(512/512)**3/(2*math.pi**2)*np.ones(np.size(pk_part_r))
    pk = (part_prop*np.sqrt(pk_part_r-shot)+(1-part_prop)*np.sqrt(pk_sl))**2
    plt.loglog(k_sl,pk,ls="-",color="green")

# strtoz = {'000': 99, '001':9,'002':4,'003':3,'004':2,'005':1,'006':0.5,'007':0.2,'008':0}
#Rounded some redshifts to nearest value...
#99, 6.5, 2.88, 1.6, 1, 0.6,0.33, 0.15, 0
strtoz = {'000': 99, '001':6.5,'002':3,'003':2,'004':1,'005':0.5,'006':0.33,'007':0.14,'008':0}

for key, value in strtoz.iteritems():
    plt.clf()
    plot_nu_powers(datadir, key,str(value),'b512p512nu0.45z99part','b512p512nu0.45z99sl')
    if value < 3:
        plot_hyb_nu_power(key, datadir, 'b512p512nu0.45z99hyb', part_prop=0.152262)
    plt.ylim(1e-11,1000)
    plt.xlabel("k /(h Mpc-1)")
    plt.ylabel("P(k) (Mpc/h)$^3$")
    save_figure(path.join(outdir,"P_nu0.45_z"+str(value)))

oonefive = {'005':1, '008':0}
for key, value in oonefive.iteritems():
    plt.clf()
    plot_nu_powers(ykdatadir, key,str(value),'part/b512p512nu0.15z49','yacine-kspace/b256p512nu0.15z99')
    plt.loglog(np.ones(np.size(np.arange(1e-12,1e4))),np.arange(1e-12,1e4),color="grey",ls="--")
    kk = key
    if value > 0.5:
        kk = '004'
    plot_hyb_nu_power(kk, datadir, 'b512p512nu0.15z99hyb', part_prop=0.0091866)
    plt.ylim(1e-11,1000)
    plt.xlabel("k /(h Mpc-1)")
    plt.ylabel("P(k) (Mpc/h)$^3$")
    save_figure(path.join(outdir,"P_nu0.15_z"+str(value)))

def make_split(nu_mass, strtoz, part_prop):
    for key, value in strtoz.iteritems():
        if value > 2:
            continue
        plt.clf()
        plt.loglog(np.ones(np.size(np.arange(1e-12,1e4))),np.arange(1e-12,1e4),color="grey",ls="--")
        plot_nu_powers(datadir, key,str(value),'b512p512nu'+nu_mass+'z99hyb','b512p512nu'+nu_mass+'z99hyb')
        plot_hyb_nu_power(key, datadir, 'b512p512nu'+nu_mass+'z99hyb', part_prop=part_prop)
        plt.ylim(1e-11,1000)
        plt.xlabel("k /(h Mpc-1)")
        plt.ylabel("P(k) (Mpc/h)$^3$")
        save_figure(path.join(outdir,"P_nu"+nu_mass+"_split_z"+str(value)))

make_split("0.15", strtoz, 0.0091866)
make_split("0.45", strtoz, 0.152262)

#High-redshift
#Total matter power comparison plots
plt.clf()
#High-redshift

def dirtofile(datadir, dirr, num_str):
    return path.join(path.join(datadir, dirr), "output/powerspec_"+num_str+".txt")

for numm in ('005','006','007','008'):
    (k, pk) = pmat.get_rel_folded_power(dirtofile(datadir, "b512p512nu0z99", numm), dirtofile(datadir, "b512p512nu0.45z99sl", numm))
    plt.semilogx(k, pk, ls="--", color="red")
    (k, pk) = pmat.get_rel_folded_power(dirtofile(datadir, "b512p512nu0z99", numm), dirtofile(datadir, "b512p512nu0.45z99part", numm))
    plt.semilogx(k, pk, ls="-.", color="black")
    (k, pk) = pmat.get_rel_folded_power(dirtofile(datadir, "b512p512nu0z99", numm), dirtofile(datadir, "b512p512nu0.45z99hyb", numm))
    plt.semilogx(k, pk, ls="-", color="green")
    plt.xlim(1e-2,10)
    plt.ylim(0.6,1)
    save_figure(path.join(outdir,"relP_tot_M045_"+numm))
    plt.clf()

