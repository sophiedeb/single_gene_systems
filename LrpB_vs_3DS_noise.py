"""LrpB_vs_3DS_noise.py"""

__author__ = "Sophie de Buyl"
__email__ = "Sophie.de.Buyl@vub.be"


from SGS import *

from matplotlib import gridspec
import matplotlib.pylab as plt
import os as ospack
from collections import OrderedDict



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

fig = plt.figure()
gs = gridspec.GridSpec(1,3)
axleg = fig.add_subplot(gs[2])
mycolors=['r', 'b', 'g', 'k']
mynonmonotonicity= 'monotonic', 'non-monotonic 1', 'non-monotonic 2', 'non-monotonic 3'


#load LrpB 3DS data
file = '/Users/sdebuyl/stoch_LrpB/'  # '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS_LrpB/'

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
counter_spikes_LrpB=0
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
tot_simulations_LrpB=len(list0)+len(list1)+len(list2)+len(list3)

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

ax = fig.add_subplot(gs[0])

myfanofactor=np.array([])

jj=1

for k in np.random.choice(len(lists[jj]),160):
    mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
    myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
    temp = mymeans_wt[vartoplot]
    if temp!=0:
        myfanofactor=np.append(myfanofactor,myvars_wt[vartoplot] / mymeans_wt[vartoplot])
        if jj==0:
            ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='0', marker="o",s=1)
        if jj==1:
            ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='1', marker="o",s=1)
        if jj==2:
            ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='2', marker="o",s=1)
        if jj==3:
            ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='3', marker="o",s=1)
        if (myvars_wt[vartoplot] / mymeans_wt[vartoplot] > 1000 or myvars_wt[vartoplot] / mymeans_wt[
            vartoplot] < .1):
            print(myvars_wt[vartoplot] / mymeans_wt[vartoplot])
        if myvars_wt[vartoplot] / mymeans_wt[vartoplot] > 2:
            counter_spikes_LrpB=counter_spikes_LrpB+1

    ax.set_xlim([0,maxnumberofmolecules])
    ax.set_ylim([0.1,1000])
    ax.set_title(mymodel+'DS LrpB')
    ax.set_ylabel('Fano factor')
    ax.set_xlabel('dimer conc.')
    ax.set_yscale('log')

print('fraction spiky sol',counter_spikes_LrpB,tot_simulations_LrpB,counter_spikes_LrpB/tot_simulations_LrpB)




########## analyse 3DS


#load LrpB 3DS data
file = '/Users/sdebuyl/stoch_' + mymodel + 'DS/'

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

myfanofactor2 = np.array([])


ax = fig.add_subplot(gs[1])

for jj in range(4):

    for k in range(len(lists[jj])):
        mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
        myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
        temp = mymeans_wt[vartoplot]
        if temp!=0:
            myfanofactor2 = np.append(myfanofactor2, myvars_wt[vartoplot] / mymeans_wt[vartoplot])
            if jj==0:
                ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='0', marker="o",s=1)
            if jj==1:
                ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='1', marker="o",s=1)
            if jj==2:
                ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='2', marker="o",s=1)
            if jj==3:
                ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='3', marker="o",s=1)
            if myvars_wt[vartoplot] / mymeans_wt[vartoplot]>1000 or  myvars_wt[vartoplot] / mymeans_wt[vartoplot]<.1:
                print( myvars_wt[vartoplot] / mymeans_wt[vartoplot])
            if myvars_wt[vartoplot] / mymeans_wt[vartoplot] > 2:
                counter_spikes_3DS = counter_spikes_3DS + 1
        ax.set_xlim([0,maxnumberofmolecules])
        ax.set_ylim([0.1,1000])
        ax.set_title(mymodel+'DS')
        ax.set_yscale('log')
        ax.set_xlabel('dimer conc.')

handles, labels = ax.get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
leg = axleg.legend(by_label.values(), by_label.keys(), ncol=1, loc=10, frameon=False, title='non-monotonicity type')
axleg.axis('off')
plt.tight_layout()


print('fraction spiky sol',counter_spikes_3DS,tot_simulations_3DS,counter_spikes_3DS/tot_simulations_3DS)


fig.savefig('noise_LrpB.pdf')

from scipy.ndimage.filters import gaussian_filter1d









#### extract data 2DS

mymodel= '2'
vartoplot=-1
oscillators = DDS


#load LrpB 3DS data
file = '/Users/sdebuyl/stoch_2DS/'  # '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS_LrpB/'

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
counter_spikes_LrpB=0
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
tot_simulations_LrpB=len(list0)+len(list1)+len(list2)+len(list3)

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



myfanofactor3=np.array([])

jj=1

for k in range(len(lists[jj])):  #np.random.choice(len(lists[jj]),160)
    mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
    myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
    temp = mymeans_wt[vartoplot]
    if temp!=0:
        myfanofactor3=np.append(myfanofactor3,myvars_wt[vartoplot] / mymeans_wt[vartoplot])



##############
##############
##############




#### extract data MDS

mymodel= 'M'
vartoplot=-1
oscillators = MDS


file = '/Users/sdebuyl/stoch_MDS/'  # '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS_LrpB/'

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
counter_spikes_LrpB=0
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
tot_simulations_LrpB=len(list0)+len(list1)+len(list2)+len(list3)

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



myfanofactor4=np.array([])

jj=1

for k in range(len(lists[jj])):  #np.random.choice(len(lists[jj]),160)
    mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
    myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
    temp = mymeans_wt[vartoplot]
    if temp!=0:
        myfanofactor4=np.append(myfanofactor4,myvars_wt[vartoplot] / mymeans_wt[vartoplot])



############## end extract MDS info
##############
##############




sorted_data4 = np.sort(myfanofactor4)
y4=np.arange(sorted_data4.size)
y4=y4.astype(float)
y4=y4/np.max(y4)

sorted_data3 = np.sort(myfanofactor3)
y3=np.arange(sorted_data3.size)
y3=y3.astype(float)
y3=y3/np.max(y3)



sorted_data = np.sort(myfanofactor)
y1=np.arange(sorted_data.size)
y1=y1.astype(float)
y1=y1/np.max(y1)
#ysmoothed = gaussian_filter1d(y1, sigma=2)
print('y1',y1)
sorted_data2 = np.sort(myfanofactor2)
y2=np.arange(sorted_data2.size)
y2=y2.astype(float)
y2=y2/np.max(y2)
fig = plt.figure()
plt.plot(sorted_data, y1/np.max(y1) ,label='LrpB')
plt.plot(sorted_data2,y2/np.max(y2) ,label='3DS')
plt.plot(sorted_data3,y3/np.max(y3),label='2DS')
plt.plot(sorted_data4,y4/np.max(y4),label='MDS')
plt.xscale('log')
#plt.xlim([0.1,10])
plt.legend()
plt.savefig('fano_factor_large_reduced.pdf')



fig = plt.figure()
plt.hist(myfanofactor4, bins='auto',alpha=.3,label='MDS')
plt.hist(myfanofactor3, bins='auto',alpha=.3,label='2DS')
plt.hist(myfanofactor2, bins='auto',alpha=.3,label='3DS')
plt.hist(myfanofactor, bins='auto',alpha=.3,label='LrpB')
plt.legend(loc='upper right')
plt.savefig('hist.pdf')

#np.savetxt('fanoLrpB.txt',myfanofactor)
#np.savetxt('fano3DS.txt',myfanofactor2)
#np.savetxt('fano2DS.txt',myfanofactor3)
#np.savetxt('fanoMDS.txt',myfanofactor4)