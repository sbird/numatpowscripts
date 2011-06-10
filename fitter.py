import numpy as np
import halofit
import glob
import os.path as path
import plot_mat_pow as pmat
import neutrinoplot as neu
import matplotlib.pyplot as plt
from scipy.optimize import leastsq,fmin
from scipy.interpolate import InterpolatedUnivariateSpline

class fitter:
        logk=np.array([])
        #Dictionary containing a tuple with (P(k), halo) objects, indexed with filename
        simpk={}
        best_fit=np.array([])
        """Args:
            linpk - File(s) containing linear theory power spectrum
            simpk - File(s) containing simulated data to fit to
            maxk - Maximum k to use in the fit
            mink - Minimum multiple of the box to use"""
        def __init__(self,infiles,maxk=6,mink=4):
            ind=np.array([])
            for s,l in infiles.iteritems():
                #Get the power
                (k,pk)=pmat.get_folded_power(s)
                #Setup self.logk if not already done
                if np.size(self.logk) ==0:
                    #Exclude the edges
                    mink=k[0]*mink
                    ind=np.where((k < maxk) * (k > mink))
                    #Rebin it to be evenly log spaced
                    self.logk=np.linspace(np.log(k[ind][0]),np.log(k[ind][-1]),np.size(k[ind])*4)
                intpk=InterpolatedUnivariateSpline(np.log(k[ind]),pk[ind])
                self.simpk[s]=(intpk(self.logk),halofit.halofit(l))
        
        """Function to find squared residuals of power spectra"""
        def halopk(self,par):
            total=0
            for s,tu in self.simpk.iteritems():
                (pk,halo)=tu
                intpk=InterpolatedUnivariateSpline(np.log(halo.k),halo.do_nonlin(par))
                halopk=intpk(self.logk)
                total+=np.sum((pk/halopk-1)**2)
            return total

        def do_min(self):
            x=fmin(self.halopk,[0,0,0])
            self.best_fit=np.array(x)
            return x

        def plot_residual(self,file,par=0):
            (pk,halo)=self.simpk[file]
            if np.size(par) < 3:
                    par=self.best_fit
            halofit=pmat.rebin(halo.do_nonlin(par),halo.k,np.exp(self.logk))
            plt.semilogx(np.exp(self.logk),pk/halofit)
            return (self.logk,pk/halofit)

#dd3=fitter.data(["/home/spb41/data3/NU_DM/PART/b150p512nu0z99","/home/spb41/data3/NU_DM/PART/b150p512nu0z99as2.0","/home/spb41/data3/NU_DM/COSMO-CHECK/b150p512nu0z99ns0.9","/home/spb41/data3/NU_DM/KSPACE/b150p512nu0z99om0.2","/home/spb41/data3/NU_DM/KSPACE/b150p512nu0z49om0.2"])

class data:
        pkfiles={} #simpk->linpk
        fit=0
        def __init__(self,simdir,pkdir="/home/spb41/cosmomc-src/cosmomc/camb/out/",maxk=6,maxz=3):
            for dir in simdir:
                files=glob.glob(path.join(dir,"powerspec_0*"))
                for f in files:
                    if neu.get_redshift(f) > maxz:
                        continue
                    self.pkfiles[f]=neu.find_linpk(f)
            self.fit=fitter(self.pkfiles,maxk=maxk)
            pars=self.fit.do_min()
            print pars



            

