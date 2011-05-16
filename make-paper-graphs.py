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

neut=neu.neutrino_power(data="/data/spb41_2/NU_DM", out="/home/spb41/Lyman-alpha/neutrinopaper/plots/")

#Kspace vs particle plots
def partvskspace(m_nu):
        if m_nu > 0.3:
                zi=str(99)
        else:
                zi=str(49)
        m_nu = str(np.around(m_nu,1))
        neut.plot_conv_diffs('/data/spb41_2/NU_DM/PART/b150p512nu'+m_nu+'z'+zi,'/data/spb41_2/NU_DM/KSPACE/b150p512nu'+m_nu+'z'+zi,ex_zz=np.array([49,9,4,3,2]))
        plt.xlim(0.05,10)
        plt.title(r'Grid vs Particle Methods, $m_\nu = '+m_nu+'$')

#512+150Mpc boxes
def halofit_sim_compare(redshift,m_nu):
        if m_nu > 0.3:
                zi=str(99)
        else:
                zi=str(49)
        m_nu=str(m_nu)
        if m_nu == "0.15":
                neut.plot_directory(['/data/spb41_2/NU_DM/KSPACE/b512p512nu'+m_nu+'z24','/data/spb41_2/NU_DM/KSPACE/b150p512nu'+m_nu+'z49'],redshifts=redshift, maxks=[1.5,])
        else:
                neut.plot_directory(['/data/spb41_2/NU_DM/KSPACE/b512p512nu'+m_nu+'z'+zi,'/data/spb41_2/NU_DM/KSPACE/b150p512nu'+m_nu+'z'+zi],redshifts=redshift, maxks=[1.5,])
        plt.xlim(0.01,10)
        zzs=re.sub(r"\.",r"_",str(np.around(redshift,2)))
        m_nus=re.sub(r"\.",r"_",m_nu)
#         plt.ylim(0.5,1)
        save_figure(path.join(neut.outdir,"nu"+m_nus+"z"+zzs))
        plt.clf()

partvskspace(0.3)
plt.ylim(-4,2)
save_figure(path.join(neut.outdir,"partvskspace0_3"))
plt.clf()
partvskspace(0.6)
plt.ylim(-6,1)
save_figure(path.join(neut.outdir,"partvskspace0_6"))
plt.clf()

for z in [0, 0.2, 0.5, 1]:
        for m in [0.15, 0.3, 0.6]:
                halofit_sim_compare(z,m)
