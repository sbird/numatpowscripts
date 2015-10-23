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
from scipy.integrate import *
from scipy.optimize import bisect,brentq
from scipy.interpolate import InterpolatedUnivariateSpline

class halofit:
    anorm = 1./(2*math.pi**2)

    def __init__(self,pkfile,omm0=0.3,omv0=0.7):
        mk1=np.loadtxt(pkfile)
        m=re.search(r"matterpow_([\d.]+).dat",pkfile)
        zz=float(m.group(1))
        m=re.search(r"om([\d.]+)_matterpow_([\d.]+).dat",pkfile)
        if m != None:
            omm0=float(m.group(1))
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

    def _neff(self,ksig,d=0):
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
        neff=self._neff(self.ksig)
        curv=self.curv(self.ksig)
        ph=self.halofit(neff,curv,self.y,par)
        pq=self.pq(self.Delta[0],self.y,neff,curv,par)
        # halo fitting formula
        return ph+pq

    def pq(self,delta,y,rn,rncur, fnu=0):
        alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
        beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ 0.3980*rn**4-0.1682*rncur
        qq=((delta+1)**beta/(delta*alpha+1))
        pq=delta*np.exp(-y/4.0-y**2/8.0)*qq
        return pq

    def halofit(self,rn,rncur,y,fnu=0):
        #Modifications from fitting halofit to the small-scale power
        gam=0.1971-0.0843*rn+0.8460*rncur
        a=1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+0.2250*rn*rn*rn*rn-0.6038*rncur
        a=10**a
        b=10**(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur)
        c=10**(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur)
        xmu=0.
        xnu=10**(5.2105+3.6902*rn)
        alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
        beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ 0.3980*rn**4-0.1682*rncur

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
        m=re.search(r"om([\d.]+)_matterpow_([\d.]+).dat",pkfile)
        if m != None:
            omm0=float(m.group(1))
        self.k=mk1[1:,0]
        self.CAMB_Pk=mk1[1:,1]/mk2[1:,1]
        self.logkmax=np.log(self.k[-1])
        self.logkmin=np.log(self.k[0])
        m=re.search(r'nu([\d.]+)',pkfile)
        self.fnu=13.4*float(m.group(1))/omm0/600.
        self.a = 1./(1+zz)
        self.om_m = self.omega_m(self.a, omm0, omv0)
        self.om_m_3 = self.omega_m(self.a, 0.3,0.7)
        self.omm0=omm0
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
        self._ph=np.empty([2,np.size(self.Delta[0])])
        self._pq=np.array(self._ph)
        self.y=np.array(self._ph)
        self.neff=np.array(self.ksig)
        self._curv=np.array(self.ksig)
        for i in [0,1]:
            ks=self.ksig[i]
            self.neff[i]=self._neff(ks,i)
            self._curv[i]=self.curv(ks,i)
            self.y[i]=self.k/ks
        self._pq[1]=self.pq(self.Delta[1],self.y[1],self.neff[1],self._curv[1],0)
        self._ph[1]=self.halofit(self.neff[1],self._curv[1],self.y[1])
        self._ph[0]=self.halofit(self.neff[0],self._curv[0],self.y[0],self.fnu)

    def do_nonlin(self,ipar):
        par=np.zeros(8)
        par[0:np.size(ipar)]=ipar
        neff=self.neff[0]
        curv=self._curv[0]
        pnl=np.array(self._ph)
        #BF from below
        par[0:4] = [  0.97737709,  47.47898431,   1.0810028,    0.39540887]
#         par[0:5] = [  0.69107114 , -0.13066274  ,47.49450963   ,1.76081586   ,0.28770219]
#         par[0:5] = [  1.38037612e+00,  2.26973755e-03,  5.79906337e+01, -1.81016587e+00, 8.70157655e-01]
#         par[0:5] = [  2.08003390e+00,   1.20133072e-03,   2.62944271e+01 , -6.48683169e+00, 1.43734127e+00]
        par[5] = ipar[0]

        self._pq[0]=self.pq(self.Delta[0],self.y[0],neff,curv,self.fnu,par)
        pnl[0]=self._pq[0]+self._ph[0]*(1+self.ph(self.a,self.ksig[0],self.fnu,neff,par))
        pnl[1]=self._pq[1]+self._ph[1]
        return pnl[0]/pnl[1]

    def pq(self,delta,y,rn,rncur, fnu,ipar=np.array([])):
        par=np.zeros(8)
        par[0:np.size(ipar)]=ipar
        alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
        beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ 0.3980*rn**4-0.1682*rncur
        beta+=fnu*(par[2]+par[3]*rn**2)
        deltaa=delta*(1+fnu*par[1]*self.k**2/(1+1.5*self.k**2))
        qq=((deltaa+1)**beta/(deltaa*alpha+1))
        pq=delta*np.exp(-y/4.0-y**2/8.0)*qq
        return pq

    def do_lin(self):
        # halo fitting formula
        return self.Delta[0]/self.Delta[1]

    def ph(self,a,ksig,fnu,rn,par):
        y=self.k
        #1.81624328e+00   5.43828928e-04(1+par[2]*(y/ksig)**0.5)
        vnl=fnu*(par[0]+par[5]*((self.omm0/0.3)**0.55-1))
        return vnl

    def halofit(self,rn,rncur,y,fnu=0,ipar=np.array([])):
        #Modifications from fitting halofit to the small-scale power
        #Fit with rescaling:

#       dd11=fitter.data(["/home/spb41/data3/NU_DM/PART/b150p512nu0z99","/home/spb41/data3/NU_DM/COSMO-CHECK/b150p512nu0z99as2.0","/home/spb41/data3/NU_DM/COSMO-CHECK/b150p512nu0z99ns0.9","/home/spb41/data3/NU_DM/COSMO-CHECK/b150p512nu0z49om0.4", "/home/spb41/data3/NU_DM/COSMO-CHECK/b150p512nu0z99h0.75"],maxz=3.1,npar=3, maxk=500,mink=6)
#Loaded files, minimising
#Optimization terminated successfully.
#         Current function value: 10.495262
#         Iterations: 7
#         Function evaluations: 257
#[ 0.31375759 -0.07979239 -0.85801047]

        # I have used a high mink because otherwise it was straining to fit the wiggles on large
        # scales. A very high maxk because otherwise it was straining to fit the jump at
        # the nyquist frequency.
        # It doesn't fit high redshift well, because that is quasilinear.
#         par=np.zeros(8)
#         par[0:np.size(ipar)]=ipar
        gam=0.1971-0.0843*rn+0.8460*rncur
        a=1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+0.2250*rn*rn*rn*rn-0.6038*rncur
        a=10**a
        b=10**(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur)
        c=10**(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur)
        xmu=0.
        xnu=10**(5.2105+3.6902*rn)
        alpha=abs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur)
        beta=2.0379-0.7354*rn+0.3157*rn**2+1.2490*rn**3+ 0.3980*rn**4-0.1682*rncur

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

#Fit, ignoring the Omega_m dependence
#This has ph = A fnu /(1+By**3)
# pq = delta*(1+fnu k**2/(1+1.5*k**2)) beta+=fnu(C+rn*D)

#dd2=fitter.data(["/home/spb41/data3/NU_DM/PART/b150p512nu0.6z99","/home/spb41/data3/NU_DM/PART/b150p512nu0.3z49","/home/spb41/data3/NU_DM/PART/b150p512nu0.15z24","/home/spb41/data3/NU_DM/PART/b512p512nu0.6z99","/home/spb41/data3/NU_DM/PART/b512p512nu0.3z49","/home/spb41/data3/NU_DM/PART/b512p512nu0.15z24"],fitter=fitter.relfit,maxz=3.1,npar=5, maxk=106,mink=6)
#Loaded files, minimising
#Optimization terminated successfully.
#         Current function value: 0.074884
#         Iterations: 15
#         Function evaluations: 840
#[  2.08003390e+00   1.20133072e-03   2.62944271e+01  -6.48683169e+00
#   1.43734127e+00]

#Fitting only for residual O_m dependence:
#dd4=fitter.data(["/home/spb41/data3/NU_DM/PART/b150p512nu0.3z49om0.25","/home/spb41/data3/NU_DM/PART/b150p512nu0.3z49"],fitter=fitter.relfit,maxz=3.1,npar=1, maxk=106,mink=6)
#Loaded files, minimising
#Optimization terminated successfully.
#         Current function value: 0.098976
#         Iterations: 2
#         Function evaluations: 37
#12.3865921838


#With Takahashi Halofit:
#dd2=fitter.data(["/home/spb/data/NU_DM/PART/b150p512nu0.6z99","/home/spb/data/NU_DM/PART/b150p512nu0.3z49","/home/spb/data/NU_DM/PART/b150p512nu0.15z24","/home/spb/data/NU_DM/PART/b512p512nu0.6z99","/home/spb/data/NU_DM/PART/b512p512nu0.3z49","/home/spb/data/NU_DM/PART/b512p512nu0.15z24"],fitter=fitter.relfit,maxz=3.1,npar=5, maxk=30,mink=0.1)
#Loaded files, minimising
#Optimization terminated successfully.
#         Current function value: 0.113300
#         Iterations: 11
#         Function evaluations: 627
#[  0.69107114  -0.13066274  47.49450963   1.76081586   0.28770219]


#dd4=fitter.data(["/home/spb/data/NU_DM/PART/b150p512nu0.3z49om0.25","/home/spb/data/NU_DM/PART/b150p512nu0.3z49","/home/spb/data/NU_DM/COSMO-CHECK/b150p512nu0.3z49om0.4"],fitter=fitter.relfit,maxz=3.1,npar=1, maxk=30,mink=0.1)
#Loaded files, minimising
#Optimization terminated successfully.
#         Current function value: 0.092392
#         Iterations: 2
#         Function evaluations: 37
#12.9598488366
#Adding the Omega_M = 0.4 simulation gives a very different result. Omega_M = 0.3, 0.4 are

#again because p[1] being negative leads to problems
#In [12]: dd2=fitter.data(["/home/spb/data/NU_DM/PART/b150p512nu0.6z99","/home/spb/data/NU_DM/PART/b150p512nu0.3z49","/home/spb/data/NU_DM/PART/b150p512nu0.15z24","/home/spb/data/NU_DM/PART/b512p512nu0.6z99","/home/spb/data/NU_DM/PART/b512p512nu0.3z49","/home/spb/data/NU_DM/PART/b512p512nu0.15z24"],fitter=fitter.relfit,maxz=3.1,npar=4, maxk=30,mink=0.1)
#Loaded files, minimising
#Optimization terminated successfully.
#         Current function value: 0.126754
#         Iterations: 7
#         Function evaluations: 306
#[  0.97737709  47.47898431   1.0810028    0.39540887]
#
# dd4=fitter.data(["/home/spb/data/NU_DM/PART/b150p512nu0.3z49om0.25","/home/spb/data/NU_DM/PART/b150p512nu0.3z49"],fitter=fitter.relfit,maxz=3.1,npar=1, maxk=30,mink=0.1)
#Loaded files, minimising
#Optimization terminated successfully.
#         Current function value: 0.091346
#         Iterations: 2
#         Function evaluations: 41
#18.0152803472
#
