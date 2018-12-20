from parametersgs import *
import itertools as it
import multiprocessing
from matplotlib import gridspec
import matplotlib.pylab as plt
import os as ospack
from collections import OrderedDict
import random as random
import matplotlib as mpl
from load_fano_from_sim import *


########################################################################################################################
# Choose analysis to be performed:
########################################################################################################################


#Plot a few time traces together with response curves and fano factors:
timetraces=False
#Plot fano factors as a function of dimer copy number
#(diff. colors for diff. types of non monotonicity) - compare the diff models
noisefig=False
#Plot figure 10 of the paper
fig10=False
#Look at multistability and high copy number (bigger than 500)
multistab=True




########################################################################################################################
# inputs to be adapted according to the model
########################################################################################################################

# choose the model ('M'=MDS,'2'=DDS, '3'=DDDS)
mymodel = '3'
# choose the model (MDS,DDS,DDDS)
oscillator = DDDS
LrpBornot = True
if fig10==True:
    LrpBornot=True

# inputs that can vary as a function of the performed analysis
maxnumberofmolecules = 10
t_f = 2000
dt = .1
tmax = t_f * dt
t = np.arange(0, t_f, dt)
# to avoid plottin too many traces on one plot:
max_number_traces_on_one_plot = 15


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


        fig = plt.figure(namefig, figsize=(15, 5))#, tight_layout=True)

        # create a subplot for the concentration of dimers as a function of time
        ax = fig.add_subplot(1, 3, 1)
        #axleg = fig.add_subplot(1, 4, 4)
        # ax.set_yscale('log')
        # 	#ax.set_ylim([0.0, 4*ssdimer])
        ax.set_xlabel('time')
        ax.set_ylabel('dimer copy number')
        # create a subplot for the response curve
        ax2 = fig.add_subplot(1, 3, 2)
        ax2.set_xlabel('monomer copy number')
        ax2.set_ylabel('activation fold')
        # create a subplot for the fano factors
        ax3 = fig.add_subplot(1,3,3)
        ax3.set_xlabel('mean dimer copy number')
        ax3.set_ylabel('fano factor')

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

            # plot a selection of time traces, and the corresponding response curve -- random.sample(max_number_traces_on_one_plot)
            #range(len(lists[jj]))
            if k==0:
                sampling=random.sample(lists[jj],max_number_traces_on_one_plot)
                #if jj==1:
                 #   [214, 120, 541, 35, 116, 431, 494, 622, 15, 1010, 815, 618, 573, 627, 403]
                tempp=lists[jj]
                print('lists jj',lists[jj])
                print('sampling',sampling)


            if np.in1d(lists[jj][k],sampling): #max_number_traces_on_one_plot):

                sample = np.arange(0,len(t),100)
                if mymeans_wt[jj][k]!=0:
                    mylabel='%.3f'%(myvars_wt[jj][k]/mymeans_wt[jj][k])

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
        #following lines where added to add legend
        #handles, labels = ax.get_legend_handles_labels()
        #by_label = OrderedDict(zip(labels, handles))
        #leg = axleg.legend(by_label.values(), by_label.keys(), ncol=1, loc=10, frameon=False, title='fano factor',
                        #fontsize=11, )
        #axleg.axis('off')
        #plt.tight_layout()

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

if noisefig==True:

     fig3 = plt.figure()
     gs = gridspec.GridSpec(2,2)

     for jj in range(4):
         ax = fig3.add_subplot(gs[jj])
         for k in range(len(lists[jj])):
             mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
             myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
             temp=mymeans_wt[vartoplot]
             if mymeans_wt[vartoplot]!=0.0:
                ax.scatter(mymeans_wt[vartoplot],myvars_wt[vartoplot]/mymeans_wt[vartoplot])
     fig3.savefig('noise.pdf')

if fig10==True:

    ########################################################################################################################
    # plot the data
    ########################################################################################################################


    # start the figure : jj labels the different types of non-monotonicity (different non-monotonicity types will be put in different figures)
    jj=1
    namefig='fig10'



    font = {'family':'Arial', 'size':10}
    mpl.rc('font', **font)
    mpl.rc('legend', handlelength=1)

    fig = plt.figure(namefig, figsize=(10,5))#, tight_layout=True)
    gs = gridspec.GridSpec(1,2, width_ratios=[5,5], height_ratios=[1])

    # create a subplot for the concentration of dimers as a function of time
    ax = fig.add_subplot(gs[0])
    ax.set_title('A', size=10, x=-0.2, y=1.1, ha='left')
    #axleg = fig.add_subplot(1, 4, 4)
    # ax.set_yscale('log')
    # 	#ax.set_ylim([0.0, 4*ssdimer])
    ax.set_xlabel('time')
    ax.set_ylabel('dimer copy number')
    # create a subplot for the histogram of fano factors
    ax2 = fig.add_subplot(gs[1])
    ax2.set_title('B', size=10, x=-0.2, y=1.1, ha='left')

    ax2.set_xlabel('Fano factor')
    ax2.set_ylabel('freq. of solutions')

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
        ts_transcient = ts#ts[:, 500:t_f]
        mymeans_wt[jj][k] = np.mean(ts_transcient[vartoplot, :])
        #print(str(lists[jj][k]), 'mean', mymeans_wt[jj][k])
        myvars_wt[jj][k] = np.var(ts_transcient[vartoplot, :])
        #print(str(lists[jj][k]), 'var', myvars_wt[jj][k])

        # plot a selection of time traces, and the corresponding response curve -- random.sample(max_number_traces_on_one_plot)
        #range(len(lists[jj]))
        if k==0:
            sampling=random.sample(lists[jj],max_number_traces_on_one_plot)
            #if jj==1:
             #   [214, 120, 541, 35, 116, 431, 494, 622, 15, 1010, 815, 618, 573, 627, 403]
            #sampling=[214, 120, 541, 35, 116, 431, 494, 622, 15, 1010, 815, 618, 573, 627, 403]
            tempp=lists[jj]
            print('lists jj',lists[jj])
            print('sampling',sampling)


        if np.in1d(lists[jj][k],sampling): #max_number_traces_on_one_plot):

            sample = np.arange(500,len(t),100)
            if mymeans_wt[jj][k]!=0:
                mylabel='%.3f'%(myvars_wt[jj][k]/mymeans_wt[jj][k])

            else:
                mylabel='nan'
            l = ax.plot(t[sample], ts[vartoplot][sample],label=mylabel,linewidth=.7)
            sample = np.arange(0, len(tt), 100)
            ax.plot(tt[sample], dde[sample,vartoplot], c=l[0].get_color(),linewidth=.7)
            # ax.set_yscale('log')
            ax.set_ylim([0, maxnumberofmolecules * 1.2])

    fanoLrpB = load_ff(mymodel='3', oscillators=DDDS, file='/Users/sdebuyl/stoch_LrpB/')
    fano3DS = load_ff(mymodel='3', oscillators=DDDS, file='/Users/sdebuyl/sgo/stoch_3DS/')
    if 1 == 9:
        fanoLrpB = load_ff(mymodel='3', oscillators=DDDS, file='/Users/sdebuyl/Dropbox/stoch_LrpB/')
        fano2DS = load_ff('2', DDS, '/Users/sdebuyl/stoch_2DS/')
        fanoMDS = load_ff('M', MDS, '/Users/sdebuyl/stoch_MDS/')
        fanoLrpB = load_ff(mymodel='3', oscillators=DDDS, file='/Users/sdebuyl/Dropbox/stoch_LrpB/')

    bins = np.linspace(0, 20, 100)
    bins = np.append(bins, np.max([np.max(fanoLrpB), np.max(fano3DS)]))
    if 1 == 1:
        # plt.hist(fanoMDS, bins=[-1.0,1,2,10,max(fanoMDS)],alpha=.3,label='MDS')
        # #plt.hist(fano2DS, bins='auto',alpha=.3,label='2DS')
        ax2.hist(fanoLrpB, bins=bins, alpha=.3, label='LrpB', normed=True)
        ax2.set_xlim([0, 20])

        ax2.hist(fano3DS, bins=bins, alpha=.3, label='3DS', normed=True)
        ax2.legend(loc='upper right')

    # print('lrpb',fanoLrpB.shape)
    print('lrpB', fanoLrpB.shape)
    # print('2ds',fano2DS.shape)
    print('3ds', fano3DS.shape)
    print('LrpB bigger than 1', len(np.where(fanoLrpB > 1.0)[0]) / len(fanoLrpB))
    print('3DS bigger than 1', len(np.where(fano3DS > 1.0)[0]) / len(fano3DS))




    fig.savefig(namefig + '.pdf')

    plt.close("all")

    # removezeros = np.where(mymeans_wt[1] != 0)[0]
    # print('number of simulations with zero mean value', len(removezeros))
    # mymeans_wt[1] = mymeans_wt[1][removezeros]
    # myvars_wt[1] = myvars_wt[1][removezeros]
    #
    # temp = [ii[0] for ii in sorted(enumerate(mymeans_wt[1]), key=lambda x: x[1])]
    #
    #
    # # mymeans_wt_ordered=mymeans_wt[1](ii)
    #
    # if LrpBornot==True:
    #     namefig2='noise_LrpB'
    # else:
    #     namefig2='noise_'+mymodel+'DS'
    #
    # fig2 = plt.figure(namefig2, figsize=(8, 5), tight_layout=True)
    #
    # plt.plot(mymeans_wt[1][temp], myvars_wt[1][temp] / mymeans_wt[1][temp], 'bo')
    # fig2.savefig(namefig2+'_noise.pdf')


if multistab==True:
    print('LrpB')
    file='/Users/sdebuyl/sgo/stoch_LrpB/'
    mylist = ospack.listdir(file)
    print('number of not too high/multi simulations ', (len(mylist)-4)/5)
    file=file +'stoch_3DS_too_high.txt'

    toohigh = np.loadtxt(file)
    print('too high', toohigh.shape)
    multi=np.loadtxt('/Users/sdebuyl/sgo/stoch_LrpB/stoch_3DS_multi_parms.txt')
    print('multi',multi.shape)
    print('')

    #stoch_3DS_multi_parms
    print('3DS')
    file='/Users/sdebuyl/sgo/stoch_3DS/'
    mylist = ospack.listdir(file)
    print('number of not too high/multi simulations ', (len(mylist)-4)/4)
    toohigh = np.loadtxt(file+'stoch_3DS_too_high.txt')
    print('too high', toohigh.shape)
    multi = np.loadtxt(file +'stoch_3DS_multi_parms.txt')
    print('multi', multi.shape)
    print('')


    #
    print('2DS')
    toohigh = np.loadtxt('/Users/sdebuyl/sgo/stoch_2DS/stoch_2DS_too_high.txt')
    print('too high', toohigh.shape)
    multi = np.loadtxt('/Users/sdebuyl/sgo/stoch_2DS/stoch_2DS_multi_parms.txt')
    print('multi', multi.shape)
    #
    print('MDS')
    toohigh = np.loadtxt('/Users/sdebuyl/stoch_MDS/stoch_MDS_too_high.txt')
    print('too high', toohigh.shape)
    multi = np.loadtxt('/Users/sdebuyl/stoch_MDS/stoch_MDS_multi_parms.txt')
    print('multi', multi.shape)







# if 2==3 :
#     if LrpBornot == True:
#         file = '/Users/sdebuyl/sgo/stoch_LrpB/'  # '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS_LrpB/'
#     else:
#         file = '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS/'
#     # LrpB_shift = 0  # this variable was introduced to be able to change label the position of numbers keeping track of different simulations
#     namefiletosavedata = 'stoch_' + mymodel + 'DS'
#     # actually import the data:
#     mylist = ospack.listdir(file)
#     print(len(mylist))
#     # create list for each type of non-monotonicity
#     list0 = list()
#     list1 = list()
#     list2 = list()
#     list3 = list()
#     LrpB_shift = 0
#     for kk in mylist:
#
#         if kk.find('stoch_' + mymodel + 'DS_parms_') != -1:
#             end = kk.find('.txt')
#             print('looking', kk[16 + LrpB_shift:17 + LrpB_shift], int(kk[18 + LrpB_shift:end]))
#             if kk[16 + LrpB_shift:17 + LrpB_shift] == '0':
#                 list0.append(int(kk[18 + LrpB_shift:end]))
#             if kk[16 + LrpB_shift:17 + LrpB_shift] == '1':
#                 list1.append(int(kk[18 + LrpB_shift:end]))
#             if kk[16 + LrpB_shift:17 + LrpB_shift] == '2':
#                 list2.append(int(kk[18 + LrpB_shift:end]))
#             if kk[16 + LrpB_shift:17 + LrpB_shift] == '3':
#                 list3.append(int(kk[18 + LrpB_shift:end]))
#
#     lists = [list0, list1, list2, list3]
#
#     print('len list 0', len(list0))
#     print('len list 1', len(list1))
#     print('len list 2', len(list2))
#     print('len list 3', len(list3))
#
#     # we will plot the fano factor as a function of the mean value of the dimer concentration
#     # we first create empty arrays to be filled with means and variances of each simulation
#     mymeans0 = np.zeros(len(list0))
#     mymeans1 = np.zeros(len(list1))
#     mymeans2 = np.zeros(len(list2))
#     mymeans3 = np.zeros(len(list3))
#
#     mymeans_wt = [mymeans0, mymeans1, mymeans2, mymeans3]
#     myvars0 = np.zeros(len(list0))
#     myvars1 = np.zeros(len(list1))
#     myvars2 = np.zeros(len(list2))
#     myvars3 = np.zeros(len(list3))
#
#     myvars_wt = [myvars0, myvars1, myvars2, myvars3]