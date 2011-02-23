import matplotlib
matplotlib.use("Agg")
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

    def plot_directory(self,dirs):
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
            print "file: "+f+" at z="+str(zz)
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
                m = re.search("nu"+m_nu,d)
                zerodir=glob.glob(d[:m.start()]+"nu0z*")
                if np.size(zerodir)>0:
                    zerodir=zerodir[0]
                else:
                    raise IOError,"No Massless neutrino simulation found for "+d
                if np.size(f2) != 0:
                    plot_rel_folded_power(path.join(zerodir,zerof),f2[0],colour=colours.pop(), ls=lss.pop())
                else:
                    print "Could not find: "+path.join(d,zerof)

            plt.title(r"P(k) change, $m_{\nu} = "+m_nu+", z="+str(zz)+"$")
            zzs=re.sub(r"\.",r"_",str(zz))
            out=path.join(outdir,"graph_z"+zzs)
            print out
            save_figure(path.join(outdir,"graph_z"+zzs))
            plt.clf()

    def get_redshift(self,file):
        f = open(file, 'r')
        a=float(f.readline())
        f.close()
        return (1./a)-1.
