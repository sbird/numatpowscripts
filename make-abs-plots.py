import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import plot_mat_pow as pmat
import numpy as np
(k_sim, delta_sim) = pmat.get_folded_power("/home/spb41/data3/NU_DM/PART/b150p512nu0z99/powerspec_008.txt")
(k_sim1, delta_sim1) = pmat.get_folded_power("/home/spb41/data3/NU_DM/PART/b150p512nu0z99/powerspec_005.txt")
plt.loglog(k_sim,delta_sim, color="green")
plt.loglog(k_sim1,delta_sim1,color="green")
(k, delta) = pmat.plot_power("/home/spb41/cosmomc-src/cosmomc/camb/out/nu0_matterpow_0.dat",0, colour="orange")
(k1, delta1) = pmat.plot_power("/home/spb41/cosmomc-src/cosmomc/camb/out/nu0_matterpow_1.dat",1, colour="orange")
plt.title("$\Lambda$CDM Matter Power spectrum at $z=0,1$")
plt.ylim(0.01,1e4)
plt.xlim(0.1,10)
plt.savefig("/home/spb41/Lyman-alpha/neutrinopaper/plots/absolute.pdf")
plt.clf()
ind = np.where(k_sim< k[-1])
relpk = pmat.rebin(delta,k, k_sim[ind])/delta_sim[ind]
plt.semilogx(k_sim[ind], relpk,color="blue")
ind = np.where(k_sim1< k1[-1])
relpk1 = pmat.rebin(delta1,k1, k_sim1[ind])/delta_sim1[ind]
plt.ylim(0.75,1.3)
plt.xlim(0.1,10)
plt.semilogx(k_sim1[ind], relpk1,ls="--", color="black")
plt.title("Halofit vs simulation at $z=0,1$")
plt.ylabel(r'Halofit / Simulation')
plt.xlabel("k /(h MPc-1)")
plt.savefig("/home/spb41/Lyman-alpha/neutrinopaper/plots/absolute-diff.pdf")
