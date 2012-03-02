import matplotlib
# matplotlib.use("TkAgg")
import re
import glob
import os
import os.path as path
import halofit as halofit_module
from plot_mat_pow import *
import matplotlib.pyplot as plt

"""Find a linear CAMB power spectrum, given a simulation.
Written for the python halofit implementation"""
def find_linpk(simpk,pkdir="/home/spb41/cosmomc-src/cosmomc/camb/out/"):
    zz=get_redshift(simpk)
    if zz >= 1 or zz <0.001:
        zz=int(zz)
    zz=np.around(zz,3)
    (m_nu, halosuf) = parse_dirname(simpk)
    linpk=path.join(pkdir,"nu"+m_nu+"-lin"+halosuf+str(zz)+".dat")
    return linpk

"""Parse the directory name for neutrino mass, particles, and starting redshift"""
def parse_dirname(dir):
    #Find the nu mass in the directory
    m=re.search(r"/b(\d+)p(\d+)nu([\d\.]+)z(\d+)",dir)
    m_nu=m.group(3)
    #Search for various other parameters we might have changed
    halosuf="_matterpow_"
    n=re.search(r"om([\d+\.]+)",dir)
    if n:
        halosuf="om"+n.group(1)+halosuf
    n=re.search(r"ns([\d+\.]+)",dir)
    if n:
        halosuf="ns"+n.group(1)+halosuf
    n=re.search(r"h([\d\.]+)",dir)
    if n:
        halosuf="h"+n.group(1)+halosuf
    n=re.search(r"as([\d+\.]+)",dir)
    if n:
        halosuf="as"+n.group(1)+halosuf
    return (m_nu, halosuf)

"""Open a file and read the first number to get its redshift"""
def get_redshift(file):
    f = open(file, 'r')
    a=float(f.readline())
    f.close()
    return (1./a)-1.

"""Find the directory corresponding to d with m_nu = 0"""
def find_zero(d):
        m=re.search(r"/b(\d+)p(\d+)nu([\d\.]+)z(\d+)",d)
        m_nu=m.group(3)
        zz=m.group(4)
        new = re.sub("nu"+m_nu,"nu0",d)
        zerodir=glob.glob(new)
        if np.size(zerodir)==0:
            print "Warning! No exact zero match found for "+d
            n=re.search("nu"+m_nu,d)
            zerodir=glob.glob(d[:n.start()]+"nu0z"+zz)
            if np.size(zerodir) == 0:
                    raise IOError,"No massless neutrino simulation found for "+d
        return zerodir[0]

"""Get the P(k) from a single file, averaging over seed values"""
def get_folded_power_with_seeds(f):
    #glob for directories with other seeds
    (d,zerof)=path.split(f)
    seeds=glob.glob(d+"see*")
    (kk, delta) = get_folded_power(path.join(d,zerof))
    if np.size(delta) == 0:
        return (kk,delta,delta)
    spk=np.array([delta,])
    #Get all power spectra
    for s in seeds:
        try:
            (kk, spka) = get_folded_power(path.join(s,zerof))
        except ValueError:
            continue
        spk=np.append(spk,[spka,],0)
    #Find mean
    total=np.shape(spk)[0]
    delta=np.sum(spk,0)/total
    disp=np.empty(0)
    #Find stddev
    if total > 1:
        disp = np.sqrt(np.sum((spk - delta)**2,0)/((total-1)*1.*total))
    return (kk, delta)

"""Get the P(k) from a single file, averaging over seed values"""
def get_pk_with_seeds(d,zerof):
    #glob for directories with other seeds
    seeds=glob.glob(d+"see*")
    (kk, relpk) = get_single_file_power(d,zerof)
    if np.size(relpk) == 0:
        return (kk,relpk,relpk)
    spk=np.array([relpk,])
    #Get all power spectra
    for s in seeds:
        try:
            (kk, spka) = get_single_file_power(s,zerof)
        except ValueError:
            continue
        spk=np.append(spk,[spka,],0)
    #Find mean
    total=np.shape(spk)[0]
    relpk=np.sum(spk,0)/total
    disp=np.empty(0)
    #Find stddev
    if total > 1:
        disp = np.sqrt(np.sum((spk - relpk)**2,0)/((total-1)*1.*total))
    return (kk, relpk,disp)

"""Get the P(k) from a single file"""
def get_single_file_power(d, zerof):
    f2 = glob.glob(path.join(d,zerof))
    zerodir=find_zero(d)
    if np.size(f2) != 0:
#         print "file: "+f2[0]
        return get_rel_folded_power(path.join(zerodir,zerof),f2[0])
    else:
        raise ValueError,"Could not find: "+path.join(d,zerof)

#Function to make all the plots in a directory.
#Arguments are: outdir, datadir (ie, where powerspec is)
class neutrino_power:
    zerodir=""
    matpowdir=""
    outdir=""
    dirs=""
    def __init__(self,out="/home/spb41/data3/NU_DM/plots",data="/home/spb41/data3/NU_DM", matpow="/home/spb41/cosmomc-src/cosmomc/camb/out/"):
        self.matpowdir=matpow
        self.outdir=out
        if data[-1] == "/":
            data = data[:-1]
        self.dirs=glob.glob(data+"/*/b*p*nu*z*")
        self.dirs = [ d for d in self.dirs if not re.search(out,d) ]
        self.dirs = [ d for d in self.dirs if not re.search("\d+nu0z\d+",d) ]
        self.dirs = [ d for d in self.dirs if not re.search("\d+z\d+seed\d*",d) ]

    def plot_halofit(self, halosuf, zz,m_nu, halofit=True):
        if halofit:
            halo=path.join(self.matpowdir,"nu0"+halosuf+str(zz)+".dat")
            if path.exists(halo):
                plt.ylabel("P(k) /(h-3 Mpc3)")
                linstyle="-."
                plt.xlabel("k /(h Mpc-1)")
                plt.title("Power spectrum change")
                (k, relpk) = get_rel_power(halo,path.join(self.matpowdir,"nu"+m_nu+halosuf+str(zz)+".dat"))
                haloff=halofit_module.halofit(halo)
                ksig = haloff.ksig
                plt.semilogx(k,relpk,color="blue")
                #Some estimate of the disp
                disp = np.log(1+k/ksig)/(1+np.log(1+k/ksig))*(float(m_nu)*0.013/0.6/0.3)
                plt.semilogx(k,relpk*(1+disp),color="green", ls="--")
                plt.semilogx(k,relpk*(1-disp),color="green", ls="--")
                linstyle="--"
            else:
                print "Could not find "+halo

    def plot_directory(self, dirs, redshifts=None, save=False, maxks=[], coloursin=None, lssin=None, halofit=True):
        if np.size(dirs) == 1:
            dirs = [dirs,]
        (m_nu, halosuf) = parse_dirname(dirs[0])
        if m_nu == '0':
            return
        files=glob.glob(dirs[0]+"/powerspec_0*")
        endpath=(re.split("/",dirs[0]))[-1]
        if save:
            outdir=path.join(self.outdir,endpath)
            if not path.exists(outdir):
                os.makedirs(outdir)
            print "Output to: "+outdir
        for f in files:
            zz=get_redshift(f)
            if redshifts != None and np.any(np.abs(redshifts-zz)/(zz+1e-7) > 0.005):
                continue
            if zz >= 1 or zz <0.001:
                zz=int(zz)
            zz=np.around(zz,3)
            zerof=path.basename(f) #Bare filename
            self.plot_halofit(halosuf,zz,m_nu,halofit=halofit) #Find halofit
            if lssin == None:
                lss=["-","-","-.",":"]
            else:
                lss=lssin
            lss.reverse()
            if coloursin == None:
                colours=["red","orange","green","black"]
            else:
                colours=coloursin
            colours.reverse()
            maxks=list(maxks)
            for d in dirs:
                (kk,relpk,disp)=get_pk_with_seeds(d,zerof)
                #plot linear theory curve from the simulation directory
                plot_rel_power(path.join(find_zero(d),"ics"+halosuf+str(zz)+".dat"),path.join(d,"ics"+halosuf+str(zz)+".dat"),colour="black", ls='--')
                if np.size(relpk) == 0:
                    continue
                if maxks != []:
                    ind=np.where(kk < maxks.pop())
                else:
                    ind=np.where(kk)
                plt.semilogx(kk[ind],relpk[ind],color=colours.pop(), ls=lss.pop())
                if kk[ind][-1] > 6:
                        print d
                        print "max dP=",(1-np.min(relpk[ind]))/(float(m_nu)*0.013/0.6/0.3),"at z=",zz,"m_nu=",m_nu

                #Find stddev
                if np.size(disp) > 0:
                    plt.semilogx(kk[ind],(relpk+disp)[ind],color="grey", ls=":")
                    plt.semilogx(kk[ind],(relpk-disp)[ind],color="grey", ls=":")

            plt.ylabel(r'$\Delta^2_\nu /\Delta^2_\mathrm{CDM}$')
            plt.xlabel("k /(h Mpc$^{-1}$)")
            plt.title(r"Change in $\Delta^2(k)$ for $M_{\nu} = "+m_nu+", z="+str(zz)+"$")
            if save:
                zzs=re.sub(r"\.",r"_",str(zz))
                plt.figure(99)
                save_figure(path.join(outdir,"graph_z"+zzs))
                plt.clf()
            elif redshifts == None:
                plt.figure()

    def plot_ics(self, dirs):
        if np.size(dirs) == 1:
            dirs = [dirs,]
        #Find the nu mass in the directory
        m=re.search(r"/b(\d+)p(\d+)nu([\d\.]+)z(\d+)",dirs[0])
        m_nu=m.group(3)
        if m_nu == '0':
            return
        f=glob.glob(path.join(dirs[0],"PK-DM-ics_*"))
        if np.size(f) != 1:
            raise IOError, "Found: "+str(np.size(f))+" IC power spectra"
        m=re.search(r"PK-DM-ics_(\d+)-(\d+)-z(\d+)-nu",f[0])
        zz=m.group(3)
        box=float(m.group(2))
        halo=path.join(self.matpowdir,"nu0_matterpow_"+zz+".dat")
        if path.exists(halo):
            plot_rel_power(halo,path.join(self.matpowdir,"nu"+m_nu+"_matterpow_"+zz+".dat"),colour="blue")
        else:
            print "Could not find "+halo
        lss=[":","-.","-"]
        colours=["grey", "green","red"]
        for d in dirs:
            zdir=self.find_zero(d)
            f = glob.glob(path.join(d,"PK-DM-ics_*"))
#             f = glob.glob(path.join(d,"powerspec_ics.txt"))
            if np.size(f) != 1:
                raise IOError, "Found: "+str(np.size(f))+" IC power spectra in "+d
            zdir=self.find_zero(d)
            f2 = glob.glob(path.join(zdir,"PK-DM-ics_*"))
#             f2 = glob.glob(path.join(zdir,"powerspec_ics.txt"))
            if np.size(f2) != 1:
                raise IOError, "Found: "+str(np.size(f2))+" IC power spectra in "+zdir
#             plot_rel_folded_power(f2[0],f[0],colour=colours.pop())
            plot_genpk_power(f2[0],f[0], box, o_nu = float(m_nu)*0.13/6, colour=colours.pop())

        plt.title(r"P(k) for ICs, $m_{\nu} = "+m_nu+", z="+str(zz)+"$")

    def plot_conv_diffs(self, dir1, dir2, ex_zz=None, lss=None):
        #Find the nu mass in the directory
        files1=map(path.basename, glob.glob(dir1+"/powerspec_0*"))
        files2=map(path.basename,glob.glob(dir2+"/powerspec_0*"))
        files = []
        #Make sure files contains only things in both directories.
        for f in files1:
            if files2.count(f) > 0:
                files.append(f)
            else:
                print "Could not find: "+path.join(dir2,f)
        if lss!=None:
            lss.reverse()
        for f in sorted(files):
            zz=get_redshift(path.join(dir1,f))
            if ex_zz != None and np.any(np.abs(np.around(ex_zz,1)-zz)/(zz+1e-7) < 0.1):
                continue
            (kk1, rpk1,disp)=get_pk_with_seeds(dir1,f)
            (kk2, rpk2,disp)=get_pk_with_seeds(dir2,f)
            rpk=100*(rebin(rpk2,kk2,kk1)-rpk1)
            if lss !=None:
                plt.semilogx(kk1,rpk,label="z="+str(np.around(zz,1)), ls=lss.pop())
            else:
                plt.semilogx(kk1,rpk,label="z="+str(np.around(zz,1)))

        plt.title(r"Change, "+dir1+"  "+dir2)
        plt.ylabel(r'$\Delta^2_\nu /\Delta^2_\mathrm{CDM}$')
        plt.xlabel("k /(h Mpc$^{-1}$)")
        plt.legend(loc=3,ncol=3, mode="expand", borderaxespad=0.)

