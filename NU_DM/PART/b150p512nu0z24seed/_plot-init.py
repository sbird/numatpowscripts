#!/usr/bin/env python

"""
Plot P(k)
"""

import numpy
import scipy.interpolate
import math
import matplotlib
matplotlib.use('agg')
from matplotlib.pyplot import *


def plot_power(box,filename_dm,filename_b,camb_filename,redshift,hub,omegab,omegam):
        scale=2.0*math.pi/box
        camb_transfer=numpy.loadtxt(camb_filename)
        sigma=2.0
        #Adjust Fourier convention.
        k=camb_transfer[:,0]*hub
        #NOW THERE IS NO h in the T anywhere.
        Pk=camb_transfer[:,1]
        #^2*2*!PI^2*2.4e-9*k*hub^3
        ylabel("$P(k) /(Mpc^3 h^{-3})$")
        xlabel("$k /(h MPc^{-1})$")
        title("Power spectrum at z="+str(redshift))
        xlim(0.01, 100)
        loglog(k/hub, Pk, linestyle="--")
        pkc=numpy.loadtxt(filename_dm)
#        final=numpy.where(pkc[:,0] < pkc[-1,0]/1.0)
#       pkc=pkc[final]
        k_gadget=(pkc[1:,0])*scale
        # The factor of 2\pi^3 corrects for the units difference between gadget and camb
        Pk_gadget_dm=(2*math.pi)**3*pkc[1:,1]/scale**3
        if ( omegab > 0 ):
                pk_b=numpy.loadtxt(filename_b)
                #               pk_b=pk_b[final]
                Pk_gadget_b=(2*math.pi)**3*pk_b[1:,1]/scale**3
        else:
                Pk_gadget_b =0

        Pk_gadget=((omegam)*Pk_gadget_dm+omegab*Pk_gadget_b)/(omegam+omegab)
        samp_err=pkc[1:,2]
        sqrt_err=numpy.array(map(math.sqrt, samp_err))
        loglog(k_gadget,Pk_gadget, color="black", linewidth="1.5")
        loglog(k_gadget,Pk_gadget*(1+sigma*(2.0/sqrt_err+1.0/samp_err)),linestyle="-.",color="black")
        loglog(k_gadget,Pk_gadget*(1-sigma*(2.0/sqrt_err+1.0/samp_err)),linestyle="-.",color="black")
        xlim(0.01,k_gadget[-1]*1.1)
        return (k_gadget,Pk_gadget)


plot_power(150,'/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/PK-DM-ics_512-150-z24-nu0.dat','/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/PK-nu-ics_512-150-z24-nu0.dat','/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/ics_matterpow_24.dat','24',0.7,0,0.25)
savefig('/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/pk-init-512-150-z24-58329.png')
clf()
plot_power(150,'/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/PK-DM-ics_512-150-z24-nu0.dat','/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/PK-nu-ics_512-150-z24-nu0.dat','/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/ics_matterpow_24.dat','24',0.7,0,0)
savefig('/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/pk-init-512-150-z24-58329-bar.png')
clf()
plot_power(150,'/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/PK-DM-ics_512-150-z24-nu0.dat','/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/PK-DM-ics_512-150-z24-nu0.dat','/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/ics_matterpow_24.dat','24',0.7,0,0.25)
savefig('/home/spb41/data3/NU_DM/PART/b150p512nu0z24seed/pk-init-512-150-z24-58329-DM.png')

