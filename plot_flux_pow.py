#!/usr/bin/env python

"""
Plot P(k)
"""

import numpy as np
import math
import scipy.interpolate
import os.path
import matplotlib
import matplotlib.pyplot as plt
import re
import glob
from save_figure import save_figure

"""Just rebins the data"""
def rebin(data, xaxis,newx):
        eps = 2e-7
#         if newx[0] < xaxis[0]-eps:
#                 raise ValueError("A value in newx is beyond the interpolation range")
        intp=scipy.interpolate.InterpolatedUnivariateSpline(np.log(xaxis),data)
        newdata=intp(np.log(newx))
        return newdata

def Hubble(zz,om,H0):
        return H0*math.sqrt(om*(1+zz)**3+(1-om))

def get_flux_power(flux,box,zz,om,H0):
        """Get the flux power from the output of the flux extractor"""
        flux_power=np.loadtxt(flux)
        bins=np.shape(flux_power)[0]
        #Units:
        #We need to divide by the box to get it into 1/Mpc units
        #and then multiply by the hubble parameter to be in 1/(km/s)
        scale=(1.0+zz)/(Hubble(zz,om,H0)*box/(H0/100.0))
        #Adjust Fourier convention.
        k=flux_power[1:,0]*scale*2.0*math.pi
        PF=flux_power[1:,1]/scale  #*((1+zz)/3.2)**3
        return (k,PF)

def get_rel_flux_power(flux1, flux2,box,zz,om,H0):
        """Get the relative flux power"""
        (k1, PF1) = get_flux_power(flux1,box,zz,om,H0)
        (k2, PF2) = get_flux_power(flux2,box,zz,om,H0)
        return (k2, PF2/PF1)

def plot_flux_power(flux_power,box,zz,om,H0,ls="-"):
        """Plot the flux power from the output of the flux extractor"""
        (k, PF) = get_flux_power(flux_power, box, zz, om, H0)
        #Adjust Fourier convention.
        plt.ylabel(r"$P_F(k) (km/s)$")
        plt.xlabel(r"$k (s/km)$")
        plt.loglog(k, PF, linestyle=ls)
        plt.xlim(k[0],k[-1])
        return

def plot_rel_flux_power(flux1,flux2,box,zz,om,H0, ls="-"):
        """Plot the flux power from the output of the flux extractor"""
        (k, rPF) = get_rel_flux_power(flux1,flux2, box, zz, om, H0)
        #Adjust Fourier convention.
        plt.xlabel(r"$k (s/km)$")
        plt.semilogx(k, rPF, linestyle=ls)
        plt.xlim(k[0],k[-1])
        return

def get_flux_power_mpc(flux,box):
        """Get the flux power from the output of the flux extractor"""
        flux_power=np.loadtxt(flux)
        bins=np.shape(flux_power)[0]
        #Units:
        #We need to divide by the box to get it into h/Mpc units
        #Adjust Fourier convention.
        k=flux_power[1:,0]*2.0*math.pi/box
        PF=flux_power[1:,1]*box
        return (k,PF)

def get_rel_flux_power_mpc(flux1, flux2,box):
        """Get the relative flux power"""
        (k1, PF1) = get_flux_power_mpc(flux1,box)
        (k2, PF2) = get_flux_power_mpc(flux2,box)
        return (k2, PF2/PF1)

def plot_rel_flux_power_mpc(flux1,flux2,box, ls="-",color="black"):
        """Plot the flux power from the output of the flux extractor"""
        (k, rPF) = get_rel_flux_power_mpc(flux1,flux2, box)
        #Adjust Fourier convention.
        plt.xlabel(r"$k$ (h/Mpc)")
        plt.ylabel(r"$\Delta P_\mathrm{F} / P_\mathrm{F}$ (%)")
        plt.semilogx(k, 100*(rPF-1), linestyle=ls,color=color)
        return

