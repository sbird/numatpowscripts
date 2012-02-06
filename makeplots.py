import matplotlib
matplotlib.use("Agg")
from power_specs import matter_pow
from plot_mat_pow import *

#Fiducial case with m_v = 0.6, Omega_nu = 0.013
#At redshift 2:
#Dashed lines are halofit
data="/home/spb41/data3/NU_NEW/"
halofit="/home/spb41/cosmomc-src/cosmomc/camb/out/"
plot_rel_power(halofit+"nu0_matterpow_2.dat",halofit+"nu0.013_matterpow_2.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_2.dat",halofit+"nu0.013-lin_matterpow_2.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512_MOREOUTPUTS/powerspec_018.txt",data+"P-Gadget3_DM_B_NU0.6_512_512_512_512/powerspec_018.txt",colour="red")
#60Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_018.txt",data+"P-Gadget3_DM_B_NU0.6_60_512_512_512/powerspec_018.txt",colour="green")
#PMGRID 1024
#plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_018.txt",data+"P-Gadget3_DM_B_NU0.6_60_512_512_1024/powerspec_018.txt",colour="green")

#Plot with my variant power spectrum estimator
#mypm=matter_pow(ob=0.05,om=0.3, base=data, H0=0.7,box=512)
#(myk,myPk)=mypm.loadpk("P-Gadget3_DM_B_512_512_512_MOREOUTPUTS/spb/PK-by-snap_018",512)
#(myk,myPknu)=mypm.loadpk("P-Gadget3_DM_B_NU0.6_512_512_512_512/spb/PK-by-snap_018",512)
#plt.semilogx(myk,myPknu/myPk*6.93875/6.57793,color="black")


#Plot the renormalised 60Mpc box
#(kk_a1,pk_a1,kk_b1,pk_b1)=loadfolded(data+"P-Gadget3_DM_B_60_512_512/powerspec_018.txt")
#(kk_a2,pk_a2,kk_b2,pk_b2)=loadfolded(data+"P-Gadget3_DM_B_NU0.6_60_512_512_512/powerspec_018.txt")
#plt.semilogx(kk_b1,rebin(pk_b2,kk_b2,kk_b1)/pk_b1/0.8*0.55,color="black")
plt.title(r"P(k) change (KSPACE), $m_{\nu} = 0.6, z=2$")
save_figure(halofit+"kspace_nu0-6_z2")
plt.clf()

#At redshift 4:
plot_rel_power(halofit+"nu0_matterpow_4.dat",halofit+"nu0.013_matterpow_4.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_4.dat",halofit+"nu0.013-lin_matterpow_4.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512_MOREOUTPUTS/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.6_512_512_512_512/powerspec_003.txt",colour="red")
#60Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.6_60_512_512_512/powerspec_003.txt",colour="green")

#plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.6_60_512_512_1024/powerspec_003.txt",colour="green")
#Plot the renormalised 60Mpc box
#(kk_a1,pk_a1,kk_b1,pk_b1)=loadfolded(data+"P-Gadget3_DM_B_60_512_512/powerspec_018.txt")
#(kk_a2,pk_a2,kk_b2,pk_b2)=loadfolded(data+"P-Gadget3_DM_B_NU0.6_60_512_512_512/powerspec_018.txt")
#semilogx(kk_b1,rebin(pk_b2,kk_b2,kk_b1)/pk_b1/0.8*0.55,color="black")

plt.title(r"P(k) change (KSPACE), $m_{\nu} = 0.6, z=4$")
save_figure(halofit+"kspace_nu0-6_z4")
plt.clf()

#At redshift 0:
#Dashed lines are halofit
data="/home/spb41/data3/NU_NEW/"
halofit="/home/spb41/cosmomc-src/cosmomc/camb/out/"
plot_rel_power(halofit+"nu0_matterpow_0.dat",halofit+"nu0.013_matterpow_0.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_0.dat",halofit+"nu0.013-lin_matterpow_0.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512_MOREOUTPUTS/powerspec_026.txt",data+"P-Gadget3_DM_B_NU0.6_512_512_512_512/powerspec_021.txt",colour="red")
plt.title(r"P(k) change (KSPACE), $m_{\nu} = 0.6, z=0$")
save_figure(halofit+"kspace_nu0-6_z0")
plt.clf()



## Same but for the particle implementation
data="/home/spb41/data3/NU_NEW/PARTS/"
halofit="/home/spb41/cosmomc-src/cosmomc/camb/out/"
plot_rel_power(halofit+"nu0_matterpow_2.dat",halofit+"nu0.013_matterpow_2.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_2.dat",halofit+"nu0.013-lin_matterpow_2.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.6_512_512_512_512/powerspec_003.txt",colour="red")
#60Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_018.txt",data+"P-Gadget3_DM_B_NU0.6_60_512_512_512/powerspec_018.txt",colour="green")

plt.title(r"P(k) change (Particles), $m_{\nu} = 0.6, z=2$")
save_figure(halofit+"parts_nu0-6_z2")
plt.clf()

#At redshift 4:
plot_rel_power(halofit+"nu0_matterpow_4.dat",halofit+"nu0.013_matterpow_4.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_4.dat",halofit+"nu0.013-lin_matterpow_4.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512/powerspec_001.txt",data+"P-Gadget3_DM_B_NU0.6_512_512_512_512/powerspec_001.txt",colour="red")
#60Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.6_60_512_512_512/powerspec_003.txt",colour="green")

plt.title(r"P(k) change (Particles), $m_{\nu} = 0.6, z=4$")
save_figure(halofit+"parts_nu0-6_z4")
plt.clf()

## Same but for the particle implementation
data="/home/spb41/data3/NU_NEW/PARTS/"
halofit="/home/spb41/cosmomc-src/cosmomc/camb/out/"
plot_rel_power(halofit+"nu0_matterpow_2.dat",halofit+"nu0.0065_matterpow_2.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_2.dat",halofit+"nu0.013-lin_matterpow_2.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.3_512_512_512_512/powerspec_003.txt",colour="red")
#60Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_018.txt",data+"P-Gadget3_DM_B_NU0.3_60_512_512_512/powerspec_018.txt",colour="green")

plt.title(r"P(k) change (Particles), $m_{\nu} = 0.3, z=2$")
save_figure(halofit+"parts_nu0-3_z2")
plt.clf()

#At redshift 4:
plot_rel_power(halofit+"nu0_matterpow_4.dat",halofit+"nu0.0065_matterpow_4.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_4.dat",halofit+"nu0.0065-lin_matterpow_4.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512/powerspec_001.txt",data+"P-Gadget3_DM_B_NU0.3_512_512_512_512/powerspec_001.txt",colour="red")
#60Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.3_60_512_512_512/powerspec_003.txt",colour="green")

plt.title(r"P(k) change (Particles), $m_{\nu} = 0.3, z=4$")
save_figure(halofit+"parts_nu0-3_z4")
plt.clf()


#m_v = 0.15, Omega_nu = 0.00325
#At redshift 2:
#Dashed lines are halofit
data="/home/spb41/data3/NU_NEW/"
halofit="/home/spb41/cosmomc-src/cosmomc/camb/out/"
plot_rel_power(halofit+"nu0_matterpow_2.dat",halofit+"nu0.0065_matterpow_2.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_2.dat",halofit+"nu0.0065-lin_matterpow_2.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512_MOREOUTPUTS/powerspec_018.txt",data+"P-Gadget3_DM_B_NU0.3_512_512_512_512/powerspec_018.txt",colour="red")
#60Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_018.txt",data+"P-Gadget3_DM_B_NU0.3_60_512_512_512/powerspec_018.txt",colour="green")
  
plt.title(r"P(k) change (KSPACE), $m_{\nu} = 0.3, z=2$")
save_figure(halofit+"kspace_nu0-15_z2")
plt.clf()

#At redshift 4:
plot_rel_power(halofit+"nu0_matterpow_4.dat",halofit+"nu0.0065_matterpow_4.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_4.dat",halofit+"nu0.0065-lin_matterpow_4.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512_MOREOUTPUTS/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.3_512_512_512_512/powerspec_003.txt",colour="red")
#60Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_60_512_512/powerspec_003.txt",data+"P-Gadget3_DM_B_NU0.3_60_512_512_512/powerspec_003.txt",colour="green")

plt.title(r"P(k) change (KSPACE), $m_{\nu} = 0.3, z=4$")
save_figure(halofit+"kspace_nu0-15_z4")
plt.clf()

#m_v = 1.2, Omega_nu = 0.026
#At redshift 2:
#Dashed lines are halofit
data="/home/spb41/data3/NU_NEW/PARTS/"
halofit="/home/spb41/cosmomc-src/cosmomc/camb/out/"
plot_rel_power(halofit+"nu0_matterpow_2.dat",halofit+"nu0.026_matterpow_2.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_2.dat",halofit+"nu0.026-lin_matterpow_2.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512/powerspec_003.txt",data+"P-Gadget3_DM_B_NU1.2_512_512_512_512/powerspec_003.txt",colour="red")
plt.title(r"P(k) change (PARTICLES), $m_{\nu} = 1.2, z=2$")
save_figure(halofit+"part_nu1-2_z2")
plt.clf()

#At redshift 4:
plot_rel_power(halofit+"nu0_matterpow_4.dat",halofit+"nu0.026_matterpow_4.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_4.dat",halofit+"nu0.026-lin_matterpow_4.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"P-Gadget3_DM_B_512_512_512/powerspec_001.txt",data+"P-Gadget3_DM_B_NU1.2_512_512_512_512/powerspec_001.txt",colour="red")
plt.title(r"P(k) change (PARTICLES), $m_{\nu} = 1.2, z=4$")
save_figure(halofit+"part_nu1-2_z4")
plt.clf()


#My KSPACE
#Same but with particles at z=49. My sims. 
#At redshift 2:
#Dashed lines are halofit
data="/home/spb41/data3/NU_DM/KSPACE/"
halofit="/home/spb41/cosmomc-src/cosmomc/camb/out/"
plot_rel_power(halofit+"nu0_matterpow_2.dat",halofit+"nu0.00325_matterpow_2.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_2.dat",halofit+"nu0.00325-lin_matterpow_2.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"b150p512nu0z99/powerspec_003.txt",data+"b150p512nu0.15z99/powerspec_003.txt",colour="red")

plt.title(r"P(k) change (KSPACE), $m_{\nu} = 0.15, z=2$")
save_figure(halofit+"kspace_nu0-15_z2z99")
plt.clf()

#At redshift 4:
plot_rel_power(halofit+"nu0_matterpow_4.dat",halofit+"nu0.00325_matterpow_4.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_4.dat",halofit+"nu0.00325-lin_matterpow_4.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"b150p512nu0z99/powerspec_001.txt",data+"b150p512nu0.15z99/powerspec_001.txt",colour="red")

plt.title(r"P(k) change (KSPACE), $m_{\nu} = 0.15, z=4$")
save_figure(halofit+"kspace_nu0-15_z4z99")
plt.clf()

#At redshift 0:
#Dashed lines are halofit
plot_rel_power(halofit+"nu0_matterpow_0.dat",halofit+"nu0.00325_matterpow_0.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_0.dat",halofit+"nu0.00325-lin_matterpow_0.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"b150p512nu0z99/powerspec_008.txt",data+"b150p512nu0.15z99/powerspec_008.txt",colour="red")
plt.title(r"P(k) change (KSPACE), $m_{\nu} = 0.15, z=0$")
save_figure(halofit+"kspace_nu0-15_z0z99")
plt.clf()

#MY PARTICLE
#Same but with particles at z=49. My sims. 
#At redshift 2:
#Dashed lines are halofit
data="/home/spb41/data3/NU_DM/PART/"
halofit="/home/spb41/cosmomc-src/cosmomc/camb/out/"
plot_rel_power(halofit+"nu0_matterpow_2.dat",halofit+"nu0.013_matterpow_2.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_2.dat",halofit+"nu0.013-lin_matterpow_2.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"b512p512nu0z49/powerspec_003.txt",data+"b512p512nu0.6z49/powerspec_003.txt",colour="red")
plot_rel_folded_power(data+"b512p512nu0z49/powerspec_003.txt",data+"b512p512nu0.6z49np1024/powerspec_003.txt",colour="orange")

plt.title(r"P(k) change (Particles), $m_{\nu} = 0.6, z=2$")
save_figure(halofit+"part_nu0-6_z2z49")
plt.clf()

#At redshift 4:
plot_rel_power(halofit+"nu0_matterpow_4.dat",halofit+"nu0.013_matterpow_4.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_4.dat",halofit+"nu0.013-lin_matterpow_4.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"b512p512nu0z49/powerspec_001.txt",data+"b512p512nu0.6z49/powerspec_001.txt",colour="red")
plot_rel_folded_power(data+"b512p512nu0z49/powerspec_001.txt",data+"b512p512nu0.6z49np1024/powerspec_001.txt",colour="orange")

plt.title(r"P(k) change (Particles), $m_{\nu} = 0.6, z=4$")
save_figure(halofit+"part_nu0-6_z4z49")
plt.clf()

#At redshift 0:
#Dashed lines are halofit
plot_rel_power(halofit+"nu0_matterpow_0.dat",halofit+"nu0.013_matterpow_0.dat")
plot_rel_power(halofit+"nu0-lin_matterpow_0.dat",halofit+"nu0.013-lin_matterpow_0.dat", colour="black")

#512Mpc box:
plot_rel_folded_power(data+"b512p512nu0z49/powerspec_008.txt",data+"b512p512nu0.6z49/powerspec_008.txt",colour="red")
plot_rel_folded_power(data+"b512p512nu0z49/powerspec_008.txt",data+"b512p512nu0.6z49np1024/powerspec_008.txt",colour="orange")
plt.title(r"P(k) change (Particles), $m_{\nu} = 0.6, z=0$")
save_figure(halofit+"part_nu0-6_z0z49")
plt.clf()


