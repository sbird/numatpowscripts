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

"""Saves the figure, automatically determining file extension"""
def save_figure(path):
        bk=matplotlib.backends.backend
        if path == "":
                return
        elif bk == 'TkAgg' or bk == 'Agg' or bk == 'GTKAgg':
                path = path+".png"
        elif bk == 'PDF' or bk == 'pdf':
                path = path+".pdf"
        elif bk == 'PS' or bk == 'ps':
                path = path+".ps"
        return plt.savefig(path)

"""Just rebins the data"""
def rebin(data, xaxis,newx):
        eps = 2e-7
#         if newx[0] < xaxis[0]-eps:
#                 raise ValueError("A value in newx is beyond the interpolation range")
        intp=scipy.interpolate.InterpolatedUnivariateSpline(np.log(xaxis),data)
        newdata=intp(np.log(newx))
        return newdata

""" Translation of my old IDL script to plot the matter power"""
def get_power(matpow_filename):
        matpow=np.loadtxt(matpow_filename)
        k=matpow[1:,0]
        Pk=matpow[1:,1]
        delta=Pk*(k/(2*math.pi))**3*4*math.pi
        #^2*2*!PI^2*2.4e-9*k*hub^3
        return(k, delta)

""" Get the neutrino power from CAMB"""
def get_nu_power(matpow_filename):
        matpow=np.loadtxt(matpow_filename)
        transfer_file=re.sub("_matterpow_","_transfer_",matpow_filename)
        
        trans=np.loadtxt(transfer_file)
        T_nu=trans[1:,5]
        T_tot=trans[1:,6]
        k=matpow[1:,0]
        Pk=matpow[1:,1]*(T_nu/T_tot)**2
        delta=Pk*(k/(2*math.pi))**3*4*math.pi
        #^2*2*!PI^2*2.4e-9*k*hub^3
        return(k, delta)

""" Plot the neutrino power from CAMB"""
def plot_nu_power(fname):
        (kk,delta)=get_nu_power(fname)
        plt.loglog(kk,delta)

""" Translation of my old IDL script to plot the matter power"""
def plot_power(matpow_filename,redshift, colour=""):
        (k,delta)=get_power(matpow_filename)
        #^2*2*!PI^2*2.4e-9*k*hub^3
        plt.ylabel(r'$\Delta$ (k)')
        plt.xlabel("k /(h Mpc-1)")
        plt.title("Power spectrum at z="+str(redshift))
        if colour == "":
            plt.loglog(k, delta, linestyle="--")
        else:
            plt.loglog(k, delta, linestyle="--",color=colour)
        return(k, delta)

"""Load a GenPk format power spectum."""
def load_genpk(path,box, o_nu = 0):
        #Load DM P(k)
        o_m = 0.3
        matpow=np.loadtxt(path)
        path_nu = re.sub("PK-DM-","PK-nu-",path)
        if glob.glob(path_nu):
                mp1a = np.loadtxt(path_nu)
                matpow_t = (mp1a*o_nu +matpow*(o_m - o_nu))/o_m
                ind = np.where(matpow_t/matpow > 2)
                matpow_t[ind] = matpow[ind]
                matpow = matpow_t
#                 raise Exception
        scale=1.0/box
        #Adjust Fourier convention.
        simk=matpow[1:,0]*scale*2.0*math.pi
        Pk=matpow[1:,1]/scale**3
        return (simk,Pk)

""" Plot the matter power as output by gen-pk"""
def plot_genpk_power(matpow1,matpow2, box,o_nu = 0, colour="blue"):
        (k, Pk1)=load_genpk(matpow1,box, o_nu)
        (k,Pk2)=load_genpk(matpow2,box, o_nu)
        #^2*2*!PI^2*2.4e-9*k*hub^3
        plt.ylabel("P(k) /(h-3 Mpc3)")
        plt.xlabel("k /(h Mpc-1)")
        plt.title("Power spectrum change")
        plt.semilogx(k, Pk2/Pk1, linestyle="-", color=colour)


""" Translation of my old IDL script to plot the matter power"""
def plot_rel_power(matpow1,matpow2, colour="blue", ls="--"):
        (k, Pk) = get_rel_power(matpow1, matpow2)
        plt.ylabel("P(k) /(h-3 Mpc3)")
        plt.xlabel("k /(h Mpc-1)")
        plt.title("Power spectrum change")
        plt.semilogx(k, Pk, linestyle=ls, color=colour)

def get_rel_power(matpow1,matpow2):
        mk1=np.loadtxt(matpow1)
        mk2=np.loadtxt(matpow2)
        k=mk1[1:,0]
        Pk1=mk1[1:,1]
        Pk2=mk2[1:,1]
        #^2*2*!PI^2*2.4e-9*k*hub^3
        return (k, Pk2/Pk1)

def get_rel_folded_power(fname1, fname2):
        #Note for some reason the small scale power is first in the file.
        (kk_a1,pk_a1,kk_b1,pk_b1)=loadfolded(fname1)
        (kk_a2,pk_a2,kk_b2,pk_b2)=loadfolded(fname2)
        relpk_a=rebin(pk_a2,kk_a2,kk_a1)/pk_a1
        relpk_b=rebin(pk_b2,kk_b2,kk_b1)/pk_b1
        #Ignore the first few bins of the b power, as they are always noisy.
        ind = np.where(kk_a1 > kk_b1[-1])
        kk_aa1 = np.ravel(kk_a1[ind])
        relpk_aa = np.ravel(relpk_a[ind])
        return (np.concatenate([kk_b1,kk_aa1]), np.concatenate([relpk_b, relpk_aa]))

def get_diff_folded_power(kk_1, relpk_1, kk_2, relpk_2):
        relpk_a=rebin(relpk_2,kk_2,kk_1)/relpk_1
        return (kk_1, relpk_a)

def plot_rel_folded_power(fname1,fname2,colour="black", ls="-"):
        (kk,relpk)=get_rel_folded_power(fname1, fname2)
        plt.semilogx(kk,relpk,color=colour, ls=ls)
        plt.ylabel(r'$\delta$ P(k)')
        plt.xlabel("k /(h Mpc-1)")
        plt.title("Power spectrum change")

def plot_folded_power(fname1):
        (kk,pk)=get_folded_power(fname1)
        plt.loglog(kk,pk)

def get_folded_power(fname1):
        (kk_a1,pk_a1,kk_b1,pk_b1)=loadfolded(fname1)
        ind = np.where(kk_a1 > kk_b1[-1])
        kk_aa1 = np.ravel(kk_a1[ind])
        pk_aa = np.ravel(pk_a1[ind])
        return (np.concatenate([kk_b1,kk_aa1]), np.concatenate([pk_b1, pk_aa]))

def calc_norm(file1, file2):
        (k1, delta1)=get_folded_power(file1)
        (k2, delta2)=get_folded_power(file2)
        max = 0.6
        min = 0.1
        ind = np.where((k1< 0.6)*(k1> 0.1))
        ind2 = np.where((k2< 0.6)*(k2> 0.1))
        diff = delta2[ind2]/rebin(delta1[ind],k1[ind],k2[ind2])
        mean=np.mean(diff)
        return mean

"""Load the neutrino power spectrum in the format output 
by the internal gadget integrator"""
def get_nu_folded_power(fname):
        f_in= np.fromfile(fname, sep=' ',count=-1)
        time=f_in[0]
        bins=f_in[1]
        data=f_in[2:(2*bins+2)].reshape(bins,2)
        scale=1000
        pk=data[:,1]
        k=data[:,0]
        delta=pk*k**3*4*math.pi
        return (k*scale, delta)

"""Plot the neutrino power spectrum in the format output 
by the internal gadget integrator"""
def plot_nu_folded_power(fname):
        (kk,delta)=get_nu_folded_power(fname)
        plt.loglog(kk,delta)

#This is a cache for files that will not change between reads
folded_filedata={}

"""Load the folded power spectrum file"""
def loadfolded(fname):
        if fname in folded_filedata and os.path.getmtime(fname) <= folded_filedata[fname][0]:
                return folded_filedata[fname][1]
        f_in= np.fromfile(fname, sep=' ',count=-1)
        #Load header
        scale=1000
        time=f_in[0]
        bins_a=f_in[1]
        mass_a=f_in[2]
        npart=f_in[3]
        #Read large scale power spectrum data
        adata=f_in[4:(10*bins_a+4)].reshape(bins_a,10)
        #Read second header
        b_off=10*bins_a+4
        time=f_in[b_off]
        bins_b=f_in[b_off+1]
        mass_b=f_in[b_off+2]
        npart=f_in[b_off+3]
        #Read small-scale data
#         print os.path.basename(fname)+' at time '+str(round((1./time-1.),2))
        bdata=f_in[(b_off+4):(10*bins_b+4+b_off)].reshape(bins_b,10)
        (kk_a, pk_a) = GetFoldedPower(adata,bins_a)
        (kk_b, pk_b) = GetFoldedPower(bdata,bins_b)
        #Ignore the sample variance dominated modes near the edge of the small-scale bins.
        if(kk_a[0] > kk_b[0]):
                ind=np.where(kk_a > 4*kk_a[0])
                kk_a=kk_a[ind]
                pk_a=pk_a[ind]
        else:
                ind=np.where(kk_b > 4*kk_b[0])
                kk_b=kk_b[ind]
                pk_b=pk_b[ind]
        folded_filedata[fname]=(os.path.getmtime(fname), (scale*kk_a, pk_a, scale*kk_b, pk_b))

        return (scale*kk_a, pk_a, scale*kk_b, pk_b)

"""Returns the dimensionless Delta parameter"""
def GetFoldedPower(adata, bins):
        #Set up variables
        # k
        K_A = adata[:,0]
        #Dimensionless power
        Delta2_A = adata[:,1]
        # Shot noise
        Shot_A = adata[:,2]
        # This is the dimensionful P[k]
        ModePow_A = adata[:,3]
        #Number of modes in a bin
        ModeCount_A = adata[:,4]
        # What is the correction? It seems to be a factor of PowerSpec_Efstathiou...
        # He divides each mode by P[k] and then after summing 
        # multiplies the whole thing by P[k] again. 
        # Is this a correction for bin centroid averaging?
        #I would normally expect this to be a window function correction, 
        #but it's too small.
        Delta2Uncorrected_A = adata[:,5]
        ModePowUncorrected_A = adata[:,6]
        #This is the Efstathiou power spectrum
        Specshape_A =  adata[:,7]
        #Total power, without dividing by the number of modes
        SumPower_A =  adata[:,8]
        # This is a volume conversion factor# 4 M_PI : [k/[2M_PI:Box]]::3
        ConvFac_A =  adata[:,9]
        MinModeCount = 20
        TargetBins = 80
        if(np.max(K_A) == 0 or np.min(K_A) ==0):
                raise Exception
        logK_A=np.log10(K_A)
        MinDlogK = (np.max(logK_A) - np.min(logK_A))/TargetBins
        istart=0
        iend=0
        k_list_A = []#np.array([])
        Pk_list = []#np.array([])
#         count_list_A =[]# np.array([])
        count=0
        targetlogK=MinDlogK+logK_A[istart]
        while(iend < bins):
                count+=ModeCount_A[iend]
                iend+=1
                if count >= MinModeCount and logK_A[iend-1] >= targetlogK:
                        pk = np.sum(SumPower_A[istart:iend])/count
                        b = np.trunc(np.sum(np.arange(istart,iend)*ModeCount_A[istart:iend])/count)
                        kk = K_A[b]
                        pk *= (ConvFac_A[b]*Specshape_A[b])
                        #This is a correction from what is done in pm_periodic. 
                        #I think it is some sort of bin weighted average.
                        #In pm_periodic he divides each mode by Specshape_A(k_mode)
                        #, and then sums them. Here he then multiplies by Specshape_A(k_bin)
                        k_list_A.append(kk)
                        Pk_list.append(pk)
#                         count_list_A.append(count)
                        istart=iend
                        targetlogK=logK_A[istart]+MinDlogK
                        count=0
                        if pk < 0:
                                raise Exception
        return (np.array(k_list_A), np.array(Pk_list))


