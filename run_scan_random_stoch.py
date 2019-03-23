"""run_scan_random_stoch.py"""

__author__ = "Sophie de Buyl"
__email__ = "Sophie.de.Buyl@vub.be"


from scan_random_stoch import *

namefiletosavedata= "./test" #"./stoch_3DS/stoch_3DS"
t_f = 2000
dt = .1
numberofcores=4
numberofparametersets=numberofcores*10
maxmolecules=500
shift=0

def main():
    # MDSn
    mainf(MDS,numberofparametersets,namefiletosavedata,t_f,dt,numberofcores,maxmolecules,shift,system=System.UNREGULATED)

    # MDS
    #mainf(MDS,numberofparametersets,namefiletosavedata,t_f,dt,numberofcores,maxmolecules,shift,system=System.RANDOM)


    #SLRPB
    #mainf(DDDS,numberofparametersets,namefiletosavedata,t_f,dt,numberofcores,maxmolecules,shift,system=System.SSLRPB)


    # 2DS
    #mainf(DDS,numberofparametersets,namefiletosavedata,t_f,dt,numberofcores,maxmolecules,shift,system=System.RANDOM)


    #3DS
    #mainf(DDDS,numberofparametersets,namefiletosavedata,t_f,dt,numberofcores,maxmolecules,shift,system=System.RANDOM)

if __name__ == "__main__":
    main()