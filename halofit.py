#!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#! The `halofit' code models the nonlinear evolution of cold matter 
#! cosmological power spectra. The full details of the way in which 
#! this is done are presented in Smith et al. (2002), MNRAS, ?, ?. 
#!
#! The code `halofit' was written by R. E. Smith & J. A. Peacock. 
#! See http://www.astro.upenn.edu/~res, 
#! Last edited 8/5/2002.
#
#! Only tested for plain LCDM models with power law initial power spectra
#
#! Adapted for F90 and CAMB, AL March 2005
#! SPB 2011: Reimplementation of halofit in python
#!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import math
import re
import numpy as np
from scipy.integrate import quad,Inf,romb
from scipy.optimize import bisect
from scipy.interpolate import InterpolatedUnivariateSpline

class halofit:
    Min_kh_nonlinear = 0.005
    om_m=0.
    om_v=0.
    CAMB_Pk = np.array([])
    intpk=0
    ksig=1e6
    k = np.array([])
    logkmax=0
    logkmin=0
    anorm = 1./(2*math.pi**2)

    def __init__(self,pkfile,omm0=0.3,omv0=0.7):
        mk1=np.loadtxt(pkfile)
        m=re.search(r"matterpow_([\d.]+).dat",pkfile)
        zz=float(m.group(1))
        self.k=mk1[1:,0]
        self.CAMB_Pk=mk1[1:,1]
        self.logkmax=np.log(self.k[-1])
        self.logkmin=np.log(self.k[0])
        a = 1./(1+zz)
        self.om_m = self.omega_m(a, omm0, omv0)  
        self.om_v = self.omega_v(a, omm0, omv0)
        # Remember => plin = k^3 * P(k) * constant
        # constant = 4*pi*V/(2*pi)^3 
        self.Delta_l = self.k**3*self.CAMB_Pk*self.anorm
#         self.intpk=InterpolatedUnivariateSpline(np.log(self.k),self.Delta_l)
        self.ksig = self._ksig()

    """Find omega_m at a given expansion factor"""
    def omega_m(self, aa,om_m0,om_v0):
        omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*aa**2.0+om_m0/aa)
        omega_m=omega_t*om_m0/(om_m0+om_v0*aa**(3.0))
        return omega_m
    
    def omega_v(self, aa,om_m0,om_v0):
        omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*aa**2.0+om_m0/aa)
        omega_v=omega_t*om_v0/(om_m0*aa**(-3)+om_v0)
        return omega_v

    # Find Delta^2(k) given k
#     def Delta(self,logk):
#         Delta=self.intpk(logk)
        #Guard against interpolation errors
#         if Pk < 0:
#             Pk = 0
#         return Delta

    #Function to pass to integrator
#     def wintegrand(self,logk,R):
#         return self.Delta(logk)*math.exp(-(np.exp(logk)*R)**2)

    def wintegrand(self,R):
        return self.Delta_l*np.exp(-(self.k*R)**2)

    def sigma2(self,logR):
        R=np.exp(logR)
        size=4*2**int(np.round(np.log2(np.size(self.k))))+1
        lspacek=np.linspace(self.logkmin,self.logkmax,size)
        intpk=InterpolatedUnivariateSpline(np.log(self.k),self.wintegrand(R))
        Delta=intpk(lspacek)
        s2 = romb(Delta, lspacek[1]-lspacek[0])
        return s2


#     def sigma2(self,logR):
#         R=np.exp(logR)
#         (s2,err) = quad(self.wintegrand,self.logkmin,self.logkmax,args=(R,),epsrel=1e-5)
#         if err > 0.01:
#                 print "Integration error: %g",err
#         return s2
    
    def sigdiff(self,logR):
        s2=self.sigma2(logR)
        return np.sqrt(s2)-1
    
    def _ksig(self):
        xlogr1=-7
        xlogr2=7
        #If no non-linear growth, return linear theory
        if self.sigdiff(xlogr1) < 0:
                return 1e6
        #Find non-linear scale k_sigma
        k_sig = 1./np.exp(bisect(self.sigdiff,xlogr1, xlogr2))
        return k_sig
    
    def neff(self,ksig):
        logR=np.log(1./ksig)
        delta = logR*0.01 #NR 5.7; cbrt(1e-6)
        return -(np.log(self.sigma2(logR+delta))-np.log(self.sigma2(logR-delta)))/(2*delta)-3

    def curv(self,ksig):
        logR=np.log(1./ksig)
        delta = logR*0.03 #NR 5.7; cbrt(1e-6)(1e-6)**0.25
        return -(np.log(self.sigma2(logR+2*delta))-2*np.log(self.sigma2(logR))+np.log(self.sigma2(logR-2*delta)))/(2*delta)**2

    """BR09 put neutrinos into the matter as well
       calculate nonlinear wavenumber (rknl), effective spectral index (rneff) and 
       curvature (rncur) of the power spectrum at the desired redshift, using method 
       described in Smith et al (2002).
       calculate nonlinear power according to halofit: pnl = pq + ph,
       where pq represents the quasi-linear (halo-halo) power and 
       where ph is represents the self-correlation halo term."""
    def do_nonlin(self,par=[0,]):
        neff=self.neff(self.ksig)
        curv=self.curv(self.ksig)
        (pnl,pq,ph) = self.halofit(neff,curv,self.ksig,self.Delta_l,par)   # halo fitting formula 
        return pnl

    """halo model nonlinear fitting formula as described in 
    Appendix C of Smith et al. (2002)"""
    def halofit(self,rn,rncur,ksig,plin,par):
        gam=par[0]+par[1]*rn+par[2]*rncur+0.86485+0.2989*rn+0.1631*rncur
        a=1.4861+1.83693*rn+1.67618*rn*rn+0.7940*rn*rn*rn+ 0.1670756*rn*rn*rn*rn-0.620695*rncur
        a=10**a
        b=10**(0.9463+0.9466*rn+0.3084*rn*rn-0.940*rncur)
        c=10**(-0.2807+0.6669*rn+0.3214*rn*rn-0.0793*rncur)
        xmu=10**(-3.54419+0.19086*rn)
        xnu=10**(0.95897+1.2857*rn)
        alpha=1.38848+0.3701*rn-0.1452*rn*rn
        beta=0.8291+0.9854*rn+0.3400*rn**2

        if abs(1-self.om_m) > 0.01: 
            f1a=self.om_m**(-0.0732)
            f2a=self.om_m**(-0.1423)
            f3a=self.om_m**(0.0725)
            f1b=self.om_m**(-0.0307)
            f2b=self.om_m**(-0.0585)
            f3b=self.om_m**(0.0743)       
            frac=self.om_v/(1.-self.om_m) 
            f1=frac*f1b + (1-frac)*f1a
            f2=frac*f2b + (1-frac)*f2a
            f3=frac*f3b + (1-frac)*f3a
        else:       
            f1=1.0
            f2=1.
            f3=1.

        y=self.k/ksig

        php=(a*y**(f1*3.))/(1+b*y**f2+(f3*c*y)**(3.-gam))   #*(1+par[0]*np.tanh((y/50)**2))
        ph=php/(1+xmu/y+xnu/y**2)
        pq=plin*((plin+1)**beta/(plin*alpha+1))*np.exp(-y/4.0-y**2/8.0)

        pnl=pq+ph
        return pnl, pq, ph



