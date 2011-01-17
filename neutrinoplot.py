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

class neutrino_power:
        zerodir=""
        matpowdir=""
        outdir=""
        dirs=""
        def __init__(self,out,data,matpow="/home/spb41/cosmomc-src/cosmomc/camb/out/"):
                self.matpowdir=matpow
                self.outdir=out
                self.dirs=glob.glob(data+"/b*p*nu*z*")
                for d in self.dirs:
                        if re.search(data+r"/b(\d+)p(\d+)nu0z(\d+)",d):
                               self.zerodir=self.dirs.pop(self.dirs.index(d))
                               break
                for d in self.dirs:
                        self.plot_directory(d)

        def plot_directory(self,dir):
                #Find the nu mass in the directory
                m=re.search(r"/b(\d+)p(\d+)nu([\d\.]+)z(\d+)",dir)
                m_nu=float(m.group(3))
                files=glob.glob(dir+"/powerspec_0*")
                for f in files:
                        zz=round(self.get_redshift(f),1)
                        if zz >= 1 or zz <0.001:
                                zz=int(zz)
                        print "file: "+f+" at z="+str(zz)
                        zerof=path.basename(f)
                        endpath=(re.split("/",path.dirname(f)))[-1]
                        halo=path.join(self.matpowdir,"nu0_matterpow_"+str(zz)+".dat")
                        if path.exists(halo):
                                plot_rel_power(halo,path.join(self.matpowdir,"nu"+str(m_nu)+"_matterpow_"+str(zz)+".dat"),colour="blue")
#                                plot_rel_power(halo,path.join(self.matpowdir,"nu"+str(m_nu)+"-lin_matterpow_"+str(zz)+".dat"), colour="black")
                        else:
                                print "Could not find "+halo
                        plot_rel_folded_power(path.join(self.zerodir,zerof),f,colour="red")
                        plt.title(r"P(k) change (KSPACE), $m_{\nu} = "+str(m_nu)+", z="+str(zz)+"$")
                        outdir=path.join(self.outdir,endpath)
                        if not path.exists(outdir):
                                os.makedirs(outdir)
                        save_figure(path.join(self.outdir,endpath+"/graph_z"+str(zz)))
                        plt.clf()

        def get_redshift(self,file):
                f = open(file, 'r')
                a=float(f.readline())
                f.close()
                return (1./a)-1.
