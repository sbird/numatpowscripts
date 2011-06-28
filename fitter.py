import numpy as np
import halofit
import glob
import os.path as path
import plot_mat_pow as pmat
import neutrinoplot as neu
import re
import matplotlib.pyplot as plt
from scipy.optimize import leastsq,fmin,fmin_powell
from scipy.interpolate import InterpolatedUnivariateSpline

class fitter:
        """Args:
            linpk - File(s) containing linear theory power spectrum
            simpk - File(s) containing simulated data to fit to
            maxk - Maximum multiple of the box to use
            mink - Minimum multiple of the box to use"""
        def __init__(self,infiles,maxk=104,mink=4):
            #Dictionary containing a tuple with (P(k), halo) objects, indexed with filename
            self.simpk={}
            for s,l in infiles.iteritems():
                #Get the power
                (k,pk)=pmat.get_folded_power(s)
                #Exclude the edges
                tmink=k[0]*mink
                tmaxk=k[0]*maxk
                ind=np.where((k < tmaxk) * (k > tmink))
                #Rebin it to be evenly log spaced
                logk=np.linspace(np.log(k[ind][0]),np.log(k[ind][-1]),np.size(k[ind])*2)
                intpk=InterpolatedUnivariateSpline(np.log(k[ind]),pk[ind])
                self.simpk[s]=(logk,intpk(logk),halofit.halofit(l))
        
        """Function to find squared residuals of power spectra"""
        def halopk(self,par):
            total=0
            for s,tu in self.simpk.iteritems():
                (logk,pk,halo)=tu
                intpk=InterpolatedUnivariateSpline(np.log(halo.k),halo.do_nonlin(par))
                halopk=intpk(logk)
                total+=np.sum((pk/halopk-1)**2)
            return total

        def do_min(self, npar=4):
            x=fmin(self.halopk,np.zeros(npar))
            self.best_fit=np.array(x)
            return x

        def plot_residual(self,file,par=0):
            (logk, pk,halo)=self.simpk[file]
            if np.size(par) < 2:
                    par=self.best_fit
            halofit=pmat.rebin(halo.do_nonlin(par),halo.k,np.exp(logk))
            plt.semilogx(np.exp(logk),pk/halofit)
        
        def plot_residual_outside(self,file,par=0):
            linpk=neu.find_linpk(file)
            (d,f)=path.split(file)
            (k,pk,disp)=neu.get_pk_with_seeds(d,f)
            ind=np.where((k <= k[0]*102) * (k >= 4*k[0]))
            #Rebin it to be evenly log spaced
            logk=np.linspace(np.log(k[ind][0]),np.log(k[ind][-1]),np.size(k[ind])*2)
            intpk=InterpolatedUnivariateSpline(np.log(k[ind]),pk[ind])
            #Rebin it to be evenly log spaced
            halo=halofit.relhalofit(linpk)
            if np.size(par) < 4:
                    par=self.best_fit
            hfit=pmat.rebin(halo.do_nonlin(par),halo.k,np.exp(logk))
            plt.semilogx(np.exp(logk),intpk(logk)/hfit)

        def plot_fit_and_sim(self,file,par=0):
            (logk,pk,halo)=self.simpk[file]
            if np.size(par) < 2:
                    par=self.best_fit
            k=np.exp(logk)
            halofit=pmat.rebin(halo.do_nonlin(par),halo.k,k)
            plt.semilogx(k,pk)
            plt.semilogx(k,halofit,ls="--")

#dd3=fitter.data(["/home/spb41/data3/NU_DM/PART/b150p512nu0z99","/home/spb41/data3/NU_DM/PART/b150p512nu0z99as2.0","/home/spb41/data3/NU_DM/COSMO-CHECK/b150p512nu0z99ns0.9","/home/spb41/data3/NU_DM/KSPACE/b150p512nu0z99om0.2","/home/spb41/data3/NU_DM/KSPACE/b150p512nu0z49om0.4"])
#Optimization terminated successfully.
#         Current function value: 30.738233
#         Iterations: 208
#         Function evaluations: 374
#[ 0.13478598 -0.30378159 -1.6199286 ]


class relfit(fitter):
        """Args:
            linpk - File(s) containing linear theory power spectrum
            simpk - File(s) containing simulated data to fit to
            maxk - Maximum multiple of the box to use
            mink - Minimum multiple of the box to use"""
        def __init__(self,infiles,maxk=102,mink=2.5):
            #Dictionary containing a tuple with (P(k), halo) objects, indexed with filename
            self.simpk={}
            for s,l in infiles.iteritems():
                #Get the power
                (d,f)=path.split(s)
                (k,relpk,disp)=neu.get_pk_with_seeds(d,f)
                #Exclude the edges
                tmink=k[0]*mink
                tmaxk=k[0]*maxk
                m=re.search(r"/b(\d+)p(\d+)nu([\d\.]+)z(\d+)",d)
                if m.group(1) == '512':
                    tmaxk=0.5
                ind=np.where((k <= tmaxk) * (k >= tmink))
                #Rebin it to be evenly log spaced
                logk=np.linspace(np.log(k[ind][0]),np.log(k[ind][-1]),np.size(k[ind])*2)
                intpk=InterpolatedUnivariateSpline(np.log(k[ind]),relpk[ind])
                self.simpk[s]=(logk,intpk(logk),halofit.relhalofit(l))
        
        def do_min(self,npar=4):
            x=fmin_powell(self.halopk,np.zeros(npar))
            self.best_fit=np.array(x)
            return x

        
class data:
        def __init__(self,simdir,fitter=fitter,pkdir="/home/spb41/cosmomc-src/cosmomc/camb/out/",maxk=104,maxz=3,npar=4,mink=4, minz=0):
            self.pkfiles={} #simpk->linpk

            for dir in simdir:
                files=glob.glob(path.join(dir,"powerspec_0*"))
                if np.size(files)==0:
                    raise ValueError, "No files found!"
                for f in files:
                    zz=neu.get_redshift(f)
                    if zz > maxz or zz < minz:
                        continue
                    self.pkfiles[f]=neu.find_linpk(f)
            self.fit=fitter(self.pkfiles,maxk=maxk,mink=mink)
            print "Loaded files, minimising"
            pars=self.fit.do_min(npar)
            print pars


            

