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
from scipy.optimize import bisect,brentq
from scipy.interpolate import InterpolatedUnivariateSpline

class halofit:
    anorm = 1./(2*math.pi**2)

    def __init__(self,pkfile,omm0=0.3,omv0=0.7):
        mk1=np.loadtxt(pkfile)
        m=re.search(r"matterpow_([\d.]+).dat",pkfile)
        zz=float(m.group(1))
        self.k=mk1[1:,0]
        self.logkmax=np.log(self.k[-1])
        self.logkmin=np.log(self.k[0])
        a = 1./(1+zz)
        self.om_m = self.omega_m(a, omm0, omv0)  
        self.om_v = self.omega_v(a, omm0, omv0)
        # Remember => plin = k^3 * P(k) * constant
        # constant = 4*pi*V/(2*pi)^3 
        self.Delta=np.empty((1,np.size(self.k)))
        self.Delta[0] = self.k**3*mk1[1:,1]*self.anorm
        self.ksig = self._ksig()
        self.y=self.k/self.ksig

    """Find omega_m at a given expansion factor"""
    def omega_m(self, aa,om_m0,om_v0):
        omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*aa**2.0+om_m0/aa)
        omega_m=omega_t*om_m0/(om_m0+om_v0*aa**(3.0))
        return omega_m
    
    def omega_v(self, aa,om_m0,om_v0):
        omega_t=1.0+(om_m0+om_v0-1.0)/(1-om_m0-om_v0+om_v0*aa**2.0+om_m0/aa)
        omega_v=omega_t*om_v0/(om_m0*aa**(-3)+om_v0)
        return omega_v

    def wintegrand(self,R,d=0):
        return self.Delta[d]*np.exp(-(self.k*R)**2)

    def sigma2(self,logR,d=0):
        R=np.exp(logR)
        size=4*2**int(np.round(np.log2(np.size(self.k))))+1
        lspacek=np.linspace(self.logkmin,self.logkmax,size)
        intpk=InterpolatedUnivariateSpline(np.log(self.k),self.wintegrand(R,d))
        Delta=intpk(lspacek)
        s2 = romb(Delta, lspacek[1]-lspacek[0])
        return s2

    def sigdiff(self,logR,d=0):
        s2=self.sigma2(logR,d)
        return np.sqrt(s2)-1
    
    def _ksig(self,d=0):
        xlogr1=-7
        xlogr2=7
        #If no non-linear growth, return linear theory
        if self.sigdiff(xlogr1) < 0:
                return 1e6
        #Find non-linear scale k_sigma
        k_sig = 1./np.exp(brentq(self.sigdiff,xlogr1, xlogr2,args=(d,)))
        return k_sig
    
    def neff(self,ksig,d=0):
        logR=np.log(1./ksig)
        delta = logR*0.01 #NR 5.7; cbrt(1e-6)
        return -(np.log(self.sigma2(logR+delta,d))-np.log(self.sigma2(logR-delta,d)))/(2*delta)-3

    def curv(self,ksig,d=0):
        logR=np.log(1./ksig)
        delta = logR*0.03 #NR 5.7; cbrt(1e-6)(1e-6)**0.25
        return -(np.log(self.sigma2(logR+2*delta,d))-2*np.log(self.sigma2(logR,d))+np.log(self.sigma2(logR-2*delta,d)))/(2*delta)**2

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
        ph=self.halofit(neff,curv,self.y,par)  
        pq=self.pq(self.Delta[0],self.y,neff,par)
        # halo fitting formula 
        return ph+pq

    def pq(self,delta,y,rn,par):
        extraalpha=par[0]*rn**7+par[1]*rn**8 #par[0]*rn**4#*rn+par[1]*rn**2
        alpha=extraalpha+1.38848+0.3701*rn-0.1452*rn*rn
        extrabeta=0#par[0]*rn**4#par[2]*rn+par[3]*rn**2#par[0]+par[1]*rn+par[2]*rn**2
        beta=extrabeta+0.8291+0.9854*rn+0.3400*rn**2
        pq=delta*((delta+1)**beta/(delta*alpha+1))*np.exp(-y/4.0-y**2/8.0)
        return pq

    """halo model nonlinear fitting formula as described in 
    Appendix C of Smith et al. (2002)"""
    def halofit(self,rn,rncur,y,par):
        extragam=0.13478598 -0.30378159*rn -1.6199286*rncur
        gam=extragam+0.86485+0.2989*rn+0.1631*rncur
        a=1.4861+1.83693*rn+1.67618*rn*rn+0.7940*rn*rn*rn+ 0.1670756*rn*rn*rn*rn-0.620695*rncur
        a=10**a
        b=10**(0.9463+0.9466*rn+0.3084*rn*rn-0.940*rncur)
        c=10**(-0.2807+0.6669*rn+0.3214*rn*rn-0.0793*rncur)
        xmu=10**(-3.54419+0.19086*rn)
        xnu=10**(0.95897+1.2857*rn)

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

        php=(a*y**(f1*3.))/(1+b*y**f2+(f3*c*y)**(3.-gam))
        ph=php/(1+xmu/y+xnu/y**2)
        return ph

class relhalofit(halofit):
    def __init__(self,pkfile,omm0=0.3,omv0=0.7):
        mk1=np.loadtxt(pkfile)
        m=re.search(r"matterpow_([\d.]+).dat",pkfile)
        #Find zero
        pklin=re.sub(r'nu([\d.]+)','nu0',pkfile)
        if pklin==pkfile:
            raise ValueError,"M_nu 0 same as M_nu != 0"
        mk2=np.loadtxt(pklin)
        zz=float(m.group(1))
        self.k=mk1[1:,0]
        self.CAMB_Pk=mk1[1:,1]/mk2[1:,1]
        self.logkmax=np.log(self.k[-1])
        self.logkmin=np.log(self.k[0])
        m=re.search(r'nu([\d.]+)',pkfile)
        self.fnu=float(m.group(1))
        self.a = 1./(1+zz)
        self.om_m = self.omega_m(self.a, omm0, omv0)  
        self.om_v = self.omega_v(self.a, omm0, omv0)
        # Remember => plin = k^3 * P(k) * constant
        # constant = 4*pi*V/(2*pi)^3 
        self.Delta=np.empty((2,np.size(self.k)))
        self.Delta[0] = self.k**3*mk1[1:,1]*self.anorm
        self.Delta[1] = self.k**3*mk2[1:,1]*self.anorm
        self.ksig = np.array([self._ksig(0),self._ksig(1)])
#         self.neff = np.array([self.neff(self.ksig[0],0),self.neff(self.ksig[1],1)])
#         self.curv = np.array([self.curv(self.ksig[0],0),self.curv(self.ksig[1],1)])
#         self.y=np.array([self.k/self.ksig[0],self.k/self.ksig[1]])
        ph=np.empty([2,np.size(self.Delta[0])])
        pq=np.array(ph)
        self.pnl=np.array(ph)
        self.rn=self.neff(self.ksig[0],0)
        for i in [0,1]:
            ks=self.ksig[i]
            neff=self.neff(ks,i)
            curv=self.curv(ks,i)
            y=self.k/ks
            ph[i]=self.halofit(neff,curv,y)  
            pq[i]=self.pq(self.Delta[i],y,neff)
            self.pnl[i]=ph[i]+pq[i]

    
    def do_nonlin(self,par):
#         y=self.k/self.ksig[0]
#         return self.Delta[0]/self.Delta[1]
        return self.pnl[0]/self.pnl[1]*(1+self.ph(self.a,self.ksig[0],self.fnu,self.rn,par))

    def pq(self,delta,y,rn):
        #Modification from fitting to the LCDM simulation, ns=0.9, as=2.0, and om=0.4
        extraalpha=3.49121094e-05*rn**11
        #OR:
#         extraalpha=-0.00411522*rn**7-0.00221331*rn**8 
        alpha=extraalpha+1.38848+0.3701*rn-0.1452*rn*rn
        beta=0.8291+0.9854*rn+0.3400*rn**2
        pq=delta*((delta+1)**beta/(delta*alpha+1))*np.exp(-y/4.0-y**2/8.0)
        return pq

    def do_lin(self):
        # halo fitting formula 
        return self.Delta[0]/self.Delta[1]
    
    def ph(self,a,ksig,fnu,rn,ipar):
        y=self.k
        par=np.zeros(7)
        da=a-0.55
        par[0:np.size(ipar)]=ipar
        ga=y-2*da
        #vnl=fnu*(y*par[0]+y**2*par[1])/(1+par[2]*(y/ksig)**2.5) #residual~0.3
        vnl=fnu*par[0]*y*(1+par[1]*y**0.5+(par[2]+par[6]*rn)*da*y+(par[3]+par[5]*rn)/y**0.8 - (par[4]+0*rn)*(ga)**2/(1+2*ga**2)) #/abs(a-0.55)**0.8)# + par[3]*ksig))# +par[4]*ksig**2))
        #Maybe add an fnu**2 term?
        return vnl

    def halofit(self,rn,rncur,y,fnu=0):
        #Modifications from fitting halofit to the small-scale power    
        extragam=0.13478598 -0.30378159*rn -1.6199286*rncur
        gam=extragam+0.86485+0.2989*rn+0.1631*rncur
        a=1.4861+1.83693*rn+1.67618*rn*rn+0.7940*rn*rn*rn+ 0.1670756*rn*rn*rn*rn-0.620695*rncur
        a=10**a
        b=10**(0.9463+0.9466*rn+0.3084*rn*rn-0.940*rncur)
        c=10**(-0.2807+0.6669*rn+0.3214*rn*rn-0.0793*rncur)
        xmu=10**(-3.54419+0.19086*rn)
        xnu=10**(0.95897+1.2857*rn)

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

        php=(a*y**(f1*3.))/(1+b*y**f2+(f3*c*y)**(3.-gam))
        ph=php/(1+xmu/y+xnu/y**2)
        return ph




