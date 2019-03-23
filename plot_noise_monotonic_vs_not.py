"""figure_noise_monotonic_vs_not.py"""

__author__ = "Sophie de Buyl"
__email__ = "Sophie.de.Buyl@vub.be"


from SGS import *

import matplotlib.pylab as plt
import os as ospack



########################################################################################################################
# inputs to be adapted according to the model
########################################################################################################################


mymodel= '3'
vartoplot=9
oscillators = DDDS
LrpBornot = False
# inputs that can vary as a function of the performed analysis
maxnumberofmolecules = 500
t_f = 2000
dt = .1
tmax = t_f * dt
t = np.arange(0, t_f, dt)

########################################################################################################################
# plot the data
########################################################################################################################


mycolors=['r', 'b', 'g', 'k']
mynonmonotonicity= 'monotonic', 'non-monotonic 1', 'non-monotonic 2', 'non-monotonic 3'

########## analyse 3DS


#load LrpB 3DS data
file = '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS/'

LrpB_shift = 0  # this variable was introduced to be able to change label the position of numbers keeping track of different simulations
namefiletosavedata = 'stoch_' + mymodel + 'DS'
# actually import the data:
mylist = ospack.listdir(file)
# create list for each type of non-monotonicity
list0 = list()
list1 = list()
list2 = list()
list3 = list()
# create a counter for spiky solution
counter_spikes_3DS=0
tot_simulations_3DS=len(list0)+len(list1)+len(list2)+len(list3)
for kk in mylist:

    if kk.find('stoch_' + mymodel + 'DS_parms_') != -1:
        end = kk.find('.txt')
        #print('looking', kk[16 + LrpB_shift:17 + LrpB_shift], int(kk[18 + LrpB_shift:end]))
        if kk[16 + LrpB_shift:17 + LrpB_shift] == '0':
            list0.append(int(kk[18 + LrpB_shift:end]))
        if kk[16 + LrpB_shift:17 + LrpB_shift] == '1':
            list1.append(int(kk[18 + LrpB_shift:end]))
        if kk[16 + LrpB_shift:17 + LrpB_shift] == '2':
            list2.append(int(kk[18 + LrpB_shift:end]))
        if kk[16 + LrpB_shift:17 + LrpB_shift] == '3':
            list3.append(int(kk[18 + LrpB_shift:end]))

lists = [list0, list1, list2, list3]

print('len list 0', len(list0))
print('len list 1', len(list1))
print('len list 2', len(list2))
print('len list 3', len(list3))

# we will plot the fano factor as a function of the mean value of the dimer concentration
# we first create empty arrays to be filled with means and variances of each simulation
mymeans0 = np.zeros(len(list0))
mymeans1 = np.zeros(len(list1))
mymeans2 = np.zeros(len(list2))
mymeans3 = np.zeros(len(list3))

mymeans_wt = [mymeans0, mymeans1, mymeans2, mymeans3]

myvars0 = np.zeros(len(list0))
myvars1 = np.zeros(len(list1))
myvars2 = np.zeros(len(list2))
myvars3 = np.zeros(len(list3))

myvars_wt = [myvars0, myvars1, myvars2, myvars3]

myfanofactor_0 = np.array([])
myfanofactor_1 = np.array([])
myfanofactor_2 = np.array([])
myfanofactor_3 = np.array([])




for jj in range(4):

    for k in range(len(lists[jj])):
        mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
        myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
        temp = mymeans_wt[vartoplot]
        if temp!=0:
            if jj==0:
                myfanofactor_0=np.append(myfanofactor_0,myvars_wt[vartoplot] / mymeans_wt[vartoplot])
            if jj==1:
                myfanofactor_1=np.append(myfanofactor_1,myvars_wt[vartoplot] / mymeans_wt[vartoplot])
            if jj==2:
                myfanofactor_2=np.append(myfanofactor_2,myvars_wt[vartoplot] / mymeans_wt[vartoplot])
            if jj==3:
                myfanofactor_3=np.append(myfanofactor_3,myvars_wt[vartoplot] / mymeans_wt[vartoplot])
            if myvars_wt[vartoplot] / mymeans_wt[vartoplot]>1000 or  myvars_wt[vartoplot] / mymeans_wt[vartoplot]<.1:
                print( myvars_wt[vartoplot] / mymeans_wt[vartoplot])
            if myvars_wt[vartoplot] / mymeans_wt[vartoplot] > 2:
                counter_spikes_3DS = counter_spikes_3DS + 1





print('fraction spiky sol',counter_spikes_3DS,tot_simulations_3DS,counter_spikes_3DS/tot_simulations_3DS)


from scipy.ndimage.filters import gaussian_filter1d









sorted_data3 = np.sort(myfanofactor_3)
y3=np.arange(sorted_data3.size)
y3=y3.astype(float)
y3=y3/np.max(y3)



sorted_data = np.sort(myfanofactor_1)
y1=np.arange(sorted_data.size)
y1=y1.astype(float)
y1=y1/np.max(y1)
#ysmoothed = gaussian_filter1d(y1, sigma=2)
print('y1',y1)

sorted_data2 = np.sort(myfanofactor_2)
y2=np.arange(sorted_data2.size)
y2=y2.astype(float)
y2=y2/np.max(y2)

sorted_data0 = np.sort(myfanofactor_0)
y0=np.arange(sorted_data0.size)
y0=y0.astype(float)
y0=y0/np.max(y0)


fig = plt.figure()
plt.plot(sorted_data, y1/np.max(y1) ,label='non monotonicity 1')
plt.plot(sorted_data2,y2/np.max(y2) ,label='non monotonicity 2')
plt.plot(sorted_data3,y3/np.max(y3),label='non monotonicity 3')
plt.plot(sorted_data0,y0/np.max(y0),label='non monotonicity 0')
plt.xscale('log')
#plt.xlim([0.1,10])
plt.legend()
plt.savefig('fano_factor_large_reduced.pdf')
