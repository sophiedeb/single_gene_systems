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
mymodel = 'M'
# choose the model (MDS,DDS,DDDS)
oscillator = MDS
LrpBornot = False
# inputs that can vary as a function of the performed analysis
maxnumberofmolecules = 500
t_f = 2000
dt = .1
tmax = t_f * dt
t = np.arange(0, t_f, dt)
# to avoid plottin too many traces on one plot:
max_number_traces_on_one_plot = 30
timetraces=False

########################################################################################################################
# inputs common to all models
########################################################################################################################
# creating an "empty" oscillator
os = oscillator()
# choose variable to plot - here dimer.
if mymodel == 'M':
    vartoplot = 4
if mymodel == '2':
    vartoplot = 5
if mymodel == '3':
    vartoplot = 9
# for dimer   MDS: 4-- DDS:5 -- DDDS:9
# physiological range of the parameters
physrange = physiologicalRange();
# parameters that are commun to all models and fixed,  (here we focussed on systems with no delay).
pars_comm = {'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'vol': 1};


########################################################################################################################
# load data, put data in lists (according to non-monotonicity type
########################################################################################################################

if LrpBornot==True:
    file = '/Users/sdebuyl/sgo/stoch_LrpB/'#'/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS_LrpB/'
else:
    file = '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS/'
#LrpB_shift = 0  # this variable was introduced to be able to change label the position of numbers keeping track of different simulations
namefiletosavedata = 'stoch_' + mymodel + 'DS'
# actually import the data:
mylist = ospack.listdir(file)
print(len(mylist))
# create list for each type of non-monotonicity
list0 = list()
list1 = list()
list2 = list()
list3 = list()
LrpB_shift=0
for kk in mylist:


    if kk.find('stoch_' + mymodel + 'DS_parms_') != -1:
        end = kk.find('.txt')
        print('looking', kk[16 + LrpB_shift:17 + LrpB_shift], int(kk[18 + LrpB_shift:end]))
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

if timetraces==True:

    ########################################################################################################################
    # plot the data
    ########################################################################################################################


    # start the figure : jj labels the different types of non-monotonicity (different non-monotonicity types will be put in different figures)
    for jj in range(4):
        if LrpBornot==True:
            namefig='timelapse_' + mymodel + 'DS_type_' + str(jj)+'_LrpB'
        else:
            namefig='timelapse_' + mymodel + 'DS_type_' + str(jj)
        fig = plt.figure(namefig, figsize=(12, 5), tight_layout=True)
        # create a subplot for the concentration of dimers as a function of time
        ax = fig.add_subplot(1, 4, 1)
        axleg = fig.add_subplot(1, 4, 4)
        # ax.set_yscale('log')
        # 	#ax.set_ylim([0.0, 4*ssdimer])
        ax.set_xlabel('time')
        ax.set_ylabel('dimer copy number')
        # create a subplot for the response curve
        ax2 = fig.add_subplot(1, 4, 2)
        # create a subplot for the fano factors
        ax3 = fig.add_subplot(1,4,3)

        # loop over all simulations for the "j" type
        for k in range(len(lists[jj])):
            myparams = np.loadtxt(file + namefiletosavedata + '_parms_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
            #mymeans = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
            #myvars = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
            ts = np.loadtxt(file + namefiletosavedata + '_ts_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
            dimts=np.shape(ts)
            pars = pars_comm
            for idx, par in enumerate(os.nameparameters):
                pars[par] = myparams[idx]
            # assign parameters "pars" to the "os" model,  and compute steady state of the model
            os.setparameters(pars)
            os.steadystate()

            # intial conditions - set in to be ss for proteins and mrna, dna starts unbound.
            myinitialcond = np.zeros(dimts[0])
            for kk in range(dimts[0]):
                myinitialcond[kk] = os.ss[os.allvars[kk]]
            #we solve the deterministic model
            tt = np.arange(0, 2000, .001)
            dde = odeint(os.eqns, myinitialcond, tt, args=(1.0, 1.0))

            # remove transient before computing mean-variance
            ts_transcient = ts #ts[:, int(t_f / 4):t_f]
            mymeans_wt[jj][k] = np.mean(ts_transcient[vartoplot, :])
            #print(str(lists[jj][k]), 'mean', mymeans_wt[jj][k])
            myvars_wt[jj][k] = np.var(ts_transcient[vartoplot, :])
            #print(str(lists[jj][k]), 'var', myvars_wt[jj][k])

            # plot a selection of time traces, and the corresponding response curve
            if (k < max_number_traces_on_one_plot):
                sample = np.arange(0,len(t),100)
                if mymeans_wt[jj][k]!=0:
                    mylabel='%.3f'%(myvars_wt[jj][k]/mymeans_wt[jj][k])
                    print(mylabel)
                else:
                    mylabel='nan'
                l = ax.plot(t[sample], ts[vartoplot][sample],label=mylabel)
                sample = np.arange(0, len(tt), 100)
                ax.plot(tt[sample], dde[sample,vartoplot], c=l[0].get_color())
                # ax.set_yscale('log')
                ax.set_ylim([0, maxnumberofmolecules * 1.2])
                xx = np.linspace(0, 1000, 1000)
                myrate = os.transcriptionrate(xx) / os.transcriptionrate(0)
                ax2.plot(xx, myrate, c=l[0].get_color())
                if mymeans_wt[jj][k] !=0:
                    ax3.scatter(mymeans_wt[jj][k],myvars_wt[jj][k]/mymeans_wt[jj][k], c=l[0].get_color())
        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        leg = axleg.legend(by_label.values(), by_label.keys(), ncol=1, loc=10, frameon=False, title='fano factor',
                           fontsize=11, )
        axleg.axis('off')
        plt.tight_layout()

        fig.savefig(namefig + '.pdf')

    plt.close("all")

    removezeros = np.where(mymeans_wt[1] != 0)[0]
    print('number of simulations with zero mean value', len(removezeros))
    mymeans_wt[1] = mymeans_wt[1][removezeros]
    myvars_wt[1] = myvars_wt[1][removezeros]

    temp = [ii[0] for ii in sorted(enumerate(mymeans_wt[1]), key=lambda x: x[1])]


    # mymeans_wt_ordered=mymeans_wt[1](ii)

    if LrpBornot==True:
        namefig2='noise_LrpB'
    else:
        namefig2='noise_'+mymodel+'DS'

    fig2 = plt.figure(namefig2, figsize=(8, 5), tight_layout=True)

    plt.plot(mymeans_wt[1][temp], myvars_wt[1][temp] / mymeans_wt[1][temp], 'bo')
    fig2.savefig(namefig2+'_noise.pdf')

else:

     fig3 = plt.figure()
     gs = gridspec.GridSpec(2,2)

     for jj in range(4):
         ax = fig3.add_subplot(gs[jj])
         for k in range(len(lists[jj])):
             mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
             myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
             temp=mymeans_wt[vartoplot]
             ax.scatter(mymeans_wt[vartoplot],myvars_wt[vartoplot]/mymeans_wt[vartoplot])
     fig3.savefig('noise.pdf')
#
#     plt.show()
# # # ax.plot(stochtemp.t,dde[:,9])
# # # ax2 = fig.add_subplot(2,1,1)
# # # ax2.plot()
# # # #fig.savefig("timelapse.pdf")
# print(fanolist)
# myindex=~np.isnan(fanolist)
# print('myindex',myindex)
# fanolist = fanolist[myindex]
# meanlist = meanlist[myindex]

# sort_indices = meanlist.argsort()
# fanolist = fanolist[sort_indices]
# meanlist = meanlist[sort_indices]

# print(meanlist)
# print(fanolist)
# fig = plt.figure('fano_'+mymodel+'DS', figsize=(8, 5), tight_layout=True)
# # #subplot cree des sous zones dans la figure
# ax = fig.add_subplot(1,1,1)
# #ax.set_yscale('log')
# # 	#ax.set_ylim([0.0, 4*ssdimer])
# ax.set_xlabel('dimer mean')
# ax.set_ylabel('dimer fano factor')
# ax.plot(meanlist,fanolist)
# plt.loglog()
# plt.show()
