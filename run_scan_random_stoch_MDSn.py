from scan_random_stoch import *


namefiletosavedata="Users/sdebuyl/stoch_MDSn/stoch_MDS"
t_f = 20
dt = .1
numberofcores=4
numberofparametersets=numberofcores*10
maxmolecules=500
shift= 0
#last variable is always False, because 2DS can not be Ss-LrpB compatible
mainf(MDS,numberofparametersets,namefiletosavedata,t_f,dt,numberofcores,maxmolecules,shift,False)

