from scan_random_stoch import *

namefiletosavedata="./stoch_2DS/stoch_2DS"
t_f = 2000
dt = .1
numberofcores=4
numberofparametersets=numberofcores*3
shift=0
maxmolecules=500

#last variable is always False, because 2DS can not be Ss-LrpB compatible
mainf(DDS,numberofparametersets,namefiletosavedata,t_f,dt,numberofcores,maxmolecules,shift,False)

