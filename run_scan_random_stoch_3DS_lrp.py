from scan_random_stoch import *

namefiletosavedata="./stoch_LrpB/stoch_3DS"
t_f = 2000
dt = .1
numberofcores=4
numberofparametersets=numberofcores*10
maxmolecules=500
shift=0
LrpBornot=True

mainf(DDDS,numberofparametersets,namefiletosavedata,t_f,dt,numberofcores,maxmolecules,shift,LrpBornot)
