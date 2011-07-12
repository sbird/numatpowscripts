import matplotlib
matplotlib.use("PDF")
from power_specs import matter_pow
import re
import glob
import os
import os.path as path
from plot_mat_pow import *
import matplotlib.pyplot as plt
import neutrinoplot as neu

neut=neu.neutrino_power(data="/data/spb41_2/NU_DM", out="/home/spb41/Lyman-alpha/neutrinopaper/plots/", matpow='/home/spb41/cosmomc-src/cosmomc/camb/out-newfit/')

#Kspace vs particle plots
#def partvskspace(m_nu):
#        if m_nu > 0.3:
#                zi="99"
#        else:
#                zi="49"
#        m_nu = str(np.around(m_nu,1))
#        neut.plot_conv_diffs('/data/spb41_2/NU_DM/PART/b150p512nu'+m_nu+'z'+zi,'/data/spb41_2/NU_DM/KSPACE/b150p512nu'+m_nu+'z'+zi,ex_zz=np.array([49,9,4,3,0.5,0.2]), lss=["-.","--","-"])
#        plt.xlim(0.05,10)
#        plt.title(r'Grid vs Particle Methods, $m_\nu = '+m_nu+'$')

#512+150Mpc boxes
def halofit_sim_compare(redshift,m_nu):
        zi = { "0.15" : "24", "0.3" : "49", "0.6" : "99"}
        m_nu=str(m_nu)
        neut.plot_directory(['/data/spb41_2/NU_DM/PART/b512p512nu'+m_nu+'z'+zi[m_nu],'/data/spb41_2/NU_DM/PART/b150p512nu'+m_nu+'z'+zi[m_nu]],redshifts=redshift, maxks=[1.5,])
        plt.xlim(0.01,10)
        plt.ylim(ymax=1.0)
        zzs=re.sub(r"\.",r"_",str(np.around(redshift,2)))
        m_nus=re.sub(r"\.",r"_",m_nu)
#         plt.ylim(0.5,1)
        save_figure(path.join(neut.outdir,"new-nu"+m_nus+"z"+zzs))
        plt.clf()

#Initial plots of grid and particle, to make the point that they are basically the same.
#neut.plot_directory(['/data/spb41_2/NU_DM/KSPACE/b150p512nu0.3z49','/data/spb41_2/NU_DM/PART/b150p512nu0.3z49'],redshifts=1, halofit=False, lssin=["-.","-"],coloursin=["red","green"])
#plt.xlim(0.02,20)
#plt.ylim(0.7,1.0)
#save_figure(path.join(neut.outdir,"shapepvsgridz1"))
#plt.clf()
#neut.plot_directory(['/data/spb41_2/NU_DM/KSPACE/b150p512nu0.3z49','/data/spb41_2/NU_DM/PART/b150p512nu0.3z49'],redshifts=0, halofit=False, lssin=["-.","-"],coloursin=["red","green"])
#plt.xlim(0.02,20)
#plt.ylim(0.7,1.0)
#save_figure(path.join(neut.outdir,"shapepvsgridz0"))
#plt.clf()

#Part vs KSPACE percentage
#partvskspace(0.3)
#plt.ylim(-5,1)
#save_figure(path.join(neut.outdir,"partvskspace0_3"))
#plt.clf()
#partvskspace(0.6)
#plt.ylim(-7,0)
#save_figure(path.join(neut.outdir,"partvskspace0_6"))
#plt.clf()

#Halofit graphs
for z in [0, 0.2, 0.5, 1,2,3]:
        for m in [0.15, 0.3, 0.6]:
                halofit_sim_compare(z,m)
