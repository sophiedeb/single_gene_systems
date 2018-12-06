from parametersgs import *
import itertools as it
import multiprocessing
from matplotlib import gridspec
import matplotlib.pylab as plt
import os as ospack
from collections import OrderedDict


########################################################################################################################
# inputs to be adapted according to the model
########################################################################################################################

# choose the model ('M'=MDS,'2'=DDS, '3'=DDDS)
mymodels = ['M','2','3']
# for dimer   MDS: 4-- DDS:5 -- DDDS:9
vartoplots=[4,5,9]
# choose the model (MDS,DDS,DDDS)
oscillators = [MDS,DDS,DDDS]
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
gs = gridspec.GridSpec(2,2)
axleg = fig.add_subplot(gs[2])
mycolors=['r', 'b', 'g', 'k']
mynonmonotonicity=['monotonic', 'non-monotonic 1', 'non-monotonic 2', 'non-monotonic 3']
for jjj in range(len(mymodels)):

    ########################################################################################################################
    # load data for the model "jjj", put data in lists (according to non-monotonicity type
    ########################################################################################################################
    mymodel = mymodels[jjj]
    vartoplot=vartoplots[jjj]
    if LrpBornot == True:
        file = '/Users/sdebuyl/sgo/stoch_LrpB/'  # '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS_LrpB/'
    else:
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
    if jjj!=2:
        ax = fig.add_subplot(gs[jjj])
    else:
        ax = fig.add_subplot(gs[3])
    for jj in range(4):

        for k in range(len(lists[jj])):
            mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
            myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
            temp = mymeans_wt[vartoplot]
            if temp!=0:
                if jj==0:
                    ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='0', marker="o",s=2)
                if jj==1:
                    ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='1', marker="o",s=2)
                if jj==2:
                    ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='2', marker="o",s=2)
                if jj==3:
                    ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=mycolors[jj], label='3', marker="o",s=2)
            ax.set_xlim([0,maxnumberofmolecules])
            ax.set_ylim([0.1,1000])
            ax.set_title(mymodel+'DS')
            ax.set_yscale('log')
            ax.set_ylabel('Fano factor')
            ax.set_xlabel('dimer conc.')
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    leg = axleg.legend(by_label.values(), by_label.keys(), ncol=1, loc=10, frameon=False, title='non-monotonicity type')
    axleg.axis('off')
    plt.tight_layout()
fig.savefig('noise.pdf')

