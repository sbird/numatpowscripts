import numpy as np
import halofit
import plot_mat_pow as pmat
from scipy.optimize import leastsq,fmin
from scipy.interpolate import InterpolatedUnivariateSpline

class fitter:
        logk=np.array([])
        simpk=np.array([])
        halo=0
        """Args:
            linpk - File containng linear theory power spectrum
            simpk - File containing simulated data to fit to
            maxk - Maximum k to use in the fit
            mink - Minimum multiple of the box to use"""
        def __init__(self,linpk,simpk,maxk=6,mink=4):
            self.halo=halofit.halofit(linpk)
            #Get the power
            (k,self.simpk)=pmat.get_folded_power(simpk)
            #Exclude the edges
            mink=k[0]*mink
            ind=np.where((k < maxk) * (k > mink))
            k=k[ind]
            #Rebin it to be evenly log spaced
            self.logk=np.linspace(np.log(k[0]),np.log(k[-1]),np.size(k)*4)
            intpk=InterpolatedUnivariateSpline(np.log(k),self.simpk[ind])
            self.simpk=intpk(self.logk)
        
        """Function to find squared residuals of power spectra"""
        def halopk(self,par):
#             return np.array([par[0], (par[1]-0.1),par[0]+par[1]-0.1])
            intpk=InterpolatedUnivariateSpline(np.log(self.halo.k),self.halo.do_nonlin(par[0],par[1]))
            halopk=intpk(self.logk)
            return np.log(halopk)-np.log(self.simpk)

        def do_min(self):
            (x,ier)=leastsq(self.halopk,[0,0])
            if ier > 4:
                print "Something wrong: ier=",ier
            return (x[0],x[1])

        
