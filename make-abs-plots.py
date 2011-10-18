import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import plot_mat_pow as pmat
import neutrinoplot as neu
import numpy as np
#Load 512 simulations
(k_sim, delta_sim) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b512p512nu0z99/powerspec_008.txt")
(k_sima, delta_sima) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b512p512nu0.6z99/powerspec_008.txt")
(k_sim1, delta_sim1) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b512p512nu0z99/powerspec_005.txt")
(k_simb, delta_simb) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b512p512nu0.6z99/powerspec_005.txt")
#Plot absolute simulations
ind=np.where(k_sim < 0.4)
plt.loglog(k_sim[ind],delta_sim[ind]/k_sim[ind], color="green")
plt.loglog(k_sima[ind],delta_sima[ind]/k_sima[ind],color="green",ls="--")
plt.loglog(k_sim1[ind],delta_sim1[ind]/k_sim1[ind], color="black")
plt.loglog(k_simb[ind],delta_simb[ind]/k_simb[ind],color="black",ls="--")
#Load simulations
(k_sim, delta_sim) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b150p512nu0z99/powerspec_008.txt")
(k_sima, delta_sima) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b150p512nu0.6z99/powerspec_008.txt")
(k_sim1, delta_sim1) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b150p512nu0z99/powerspec_005.txt")
(k_simb, delta_simb) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b150p512nu0.6z99/powerspec_005.txt")
#Load linear theory
(k, delta_lin) = pmat.get_power("/home/spb41/cosmomc-src/cosmomc/camb/out/nu0-lin_matterpow_0.dat")
(k1, delta1_lin) = pmat.get_power("/home/spb41/cosmomc-src/cosmomc/camb/out/nu0-lin_matterpow_1.dat")
(ka, delta_lina) = pmat.get_power("/home/spb41/cosmomc-src/cosmomc/camb/out/nu0.6-lin_matterpow_0.dat")
(k1a, delta1_lina) = pmat.get_power("/home/spb41/cosmomc-src/cosmomc/camb/out/nu0.6-lin_matterpow_1.dat")
plt.loglog(k,delta_lin/k, color="green",ls="-.")
plt.loglog(ka,delta_lina/k,color="green",ls=":")
plt.loglog(k1,delta1_lin/k, color="black",ls="-.")
plt.loglog(k1a,delta1_lina/k,color="black",ls=":")
ind=np.where(k_sim > k_sim[0]*4)
plt.loglog(k_sim[ind],delta_sim[ind]/k_sim[ind], color="green")
plt.loglog(k_sima[ind],delta_sima[ind]/k_sima[ind],color="green",ls="--")
plt.loglog(k_sim1[ind],delta_sim1[ind]/k_sim1[ind], color="black")
plt.loglog(k_simb[ind],delta_simb[ind]/k_simb[ind],color="black",ls="--")
#Setup plot titles
plt.title("Matter Power Spectra at $z=0,1$")
plt.ylabel(r'$\Delta^2$(k)/k /(Mpc h-1)')
plt.xlabel("k /(h Mpc-1)")
plt.ylim(0.5,1e2)
plt.xlim(0.05,10)
plt.savefig("/home/spb41/Lyman-alpha/neutrinopaper/plots/absolute.pdf")
plt.clf()
#Get 512Mpc boxes
(k_big, delta_big) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b512p512nu0z99/powerspec_008.txt")
(k_big1, delta_big1) = neu.get_folded_power_with_seeds("/home/spb41/data3/NU_DM/PART/b512p512nu0z99/powerspec_005.txt")
#Load halofit
(k, delta) = pmat.get_power("/home/spb41/cosmomc-src/cosmomc/camb/out/nu0_matterpow_0.dat")
# plt.loglog(k,delta/k,ls="--", color="orange")
(k1, delta1) = pmat.get_power("/home/spb41/cosmomc-src/cosmomc/camb/out/nu0_matterpow_1.dat")
#Plot relative simulations
#Rebin
ind = np.where((k_sim< k[-1])) #*(k_sim > k_sim[0]*6) 
relpk = pmat.rebin(delta,k, k_sim[ind])/delta_sim[ind]
plt.semilogx(k_sim[ind], relpk,color="blue")
ind = np.where(k_big< 0.4)
relpk_big = pmat.rebin(delta,k, k_big[ind])/delta_big[ind]
plt.semilogx(k_big[ind], relpk_big,color="blue")
ind = np.where((k_sim1< k1[-1]))  #*(k_sim1 > k_sim1[0]*6)
relpk1 = pmat.rebin(delta1,k1, k_sim1[ind])/delta_sim1[ind]
#Setup plot titles
plt.semilogx(k_sim1[ind], relpk1,ls="--", color="black")
ind = np.where(k_big1< 0.4)
relpk_big1 = pmat.rebin(delta1,k, k_big1[ind])/delta_big1[ind]
plt.semilogx(k_big1[ind], relpk_big1,ls="--", color="black")
plt.ylim(0.75,1.3)
plt.xlim(0.01,10)
plt.title("Halofit vs simulation at $z=0,1$")
plt.ylabel(r'Halofit / Simulation')
plt.xlabel("k /(h Mpc-1)")
plt.savefig("/home/spb41/Lyman-alpha/neutrinopaper/plots/absolute-diff.pdf")
plt.clf()

#Make absolute comparison plot
ind = np.where((k_sim< k[-1])) #*(k_sim > k_sim[0]*6) 
relpk = pmat.rebin(delta,k, k_sim[ind])/delta_sim[ind]
plt.loglog(k, delta,ls="--", color="blue")
plt.loglog(k_sim[ind], delta_sim[ind],ls="-", color="black")
plt.ylim(1e-3,1e3)
plt.xlim(0.01,10)
plt.title("Halofit vs simulation at $z=0$")
plt.ylabel(r'$\Delta^2$(k)/k /(Mpc h-1)')
plt.xlabel("k /(h Mpc-1)")
plt.savefig("/home/spb41/Lyman-alpha/neutrinopaper/plots/absolute-compare.pdf")
