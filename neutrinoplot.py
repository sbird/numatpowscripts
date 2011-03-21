import matplotlib
matplotlib.use("TkAgg")
from power_specs import matter_pow
import re
import glob
import os
import os.path as path
from plot_mat_pow import *
import matplotlib.pyplot as plt


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

    def plot_directory(self,dirs, save=False):
        if np.size(dirs) == 1:
            dirs = [dirs,]
        #Find the nu mass in the directory
        m=re.search(r"/b(\d+)p(\d+)nu([\d\.]+)z(\d+)",dirs[0])
        m_nu=m.group(3)
        if m_nu == '0':
            return
        files=glob.glob(dirs[0]+"/powerspec_0*")
        endpath=(re.split("/",dirs[0]))[-1]
        outdir=path.join(self.outdir,endpath)
        if not path.exists(outdir):
            os.makedirs(outdir)
        print "Output to: "+outdir
        for f in files:
            zz=round(self.get_redshift(f),1)
            if zz >= 1 or zz <0.001:
                zz=int(zz)
            zerof=path.basename(f)
            halo=path.join(self.matpowdir,"nu0_matterpow_"+str(zz)+".dat")
            if path.exists(halo):
                plot_rel_power(halo,path.join(self.matpowdir,"nu"+m_nu+"_matterpow_"+str(zz)+".dat"),colour="blue")
#                plot_rel_power(halo,path.join(self.matpowdir,"nu"+str(m_nu)+"-lin_matterpow_"+str(zz)+".dat"), colour="black")
            else:
                print "Could not find "+halo
            lss=[":","-.","-"]
            colours=["grey", "green","red"]
            for d in dirs:
                f2 = glob.glob(path.join(d,zerof))
                zerodir=self.find_zero(d)
                if np.size(f2) != 0:
                    print "file: "+f2[0]+" at z="+str(zz)
                    plot_rel_folded_power(path.join(zerodir,zerof),f2[0],colour=colours.pop(), ls=lss.pop())
                else:
                    print "Could not find: "+path.join(d,zerof)

            plt.title(r"P(k) change, $m_{\nu} = "+m_nu+", z="+str(zz)+"$")
            zzs=re.sub(r"\.",r"_",str(zz))
            if save:
                plt.figure(99)
                save_figure(path.join(outdir,"graph_z"+zzs))
                plt.clf()
            else:
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

    def find_zero(self, d):
            m=re.search(r"/b(\d+)p(\d+)nu([\d\.]+)z(\d+)",d)
            m_nu=m.group(3)
            zz=m.group(4)
            new = re.sub("nu"+m_nu,"nu0",d)
            zerodir=glob.glob(new)
            if np.size(zerodir)==0:
                n=re.search("nu"+m_nu,d)
                zerodir=glob.glob(d[:n.start()]+"nu0z"+zz)
                if np.size(zerodir) == 0:
                        raise IOError,"No massless neutrino simulation found for "+d
            return zerodir[0]
    
    def plot_conv_diffs(self, dir1, dir2):
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
        zdir1=self.find_zero(dir1)
        zdir2=self.find_zero(dir2)
        line=np.array([])
        legname=np.array([])
        for f in files:
            (kk1, rpk1)=get_rel_folded_power(path.join(zdir1,f),path.join(dir1,f))
            (kk2, rpk2)=get_rel_folded_power(path.join(zdir2,f),path.join(dir2,f))
            (kk, rpk)=get_diff_folded_power(kk1,rpk1,kk2,rpk2)
            line=np.append(line,plt.semilogx(kk,rpk))
            legname=np.append(legname,f)

        plt.title(r"Change, "+dir1+"  "+dir2)
        plt.ylabel(r'Percentage change')
        plt.xlabel("k /(h MPc-1)")
        plt.legend(line,legname)

    def get_redshift(self,file):
        f = open(file, 'r')
        a=float(f.readline())
        f.close()
        return (1./a)-1.
