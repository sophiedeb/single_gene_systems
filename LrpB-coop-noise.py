"""LrpB-coop-noise.py"""

__author__ = "Sophie de Buyl"
__email__ = "Sophie.de.Buyl@vub.be"


from SGS import *
from matplotlib import gridspec
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

fig = plt.figure()
gs = gridspec.GridSpec(1,2,width_ratios=[8,1])
axleg = fig.add_subplot(gs[1])
mycolors=['r', 'b', 'g', 'k']
mynonmonotonicity= 'monotonic', 'non-monotonic 1', 'non-monotonic 2', 'non-monotonic 3'


#load LrpB 3DS data
file = '/Users/sdebuyl/sgo-py/stoch_3DS/'  # '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS_LrpB/'

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

ax = fig.add_subplot(gs[0])

myparamsn = oscillators().nameparameters

#  3DS np.array(['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0',
#                                        'f1', 'f2', 'f3','f12', 'f13','f23','f123',
#                                        'kb1','kb2','kb3','ku1', 'ku2', 'ku3',
#                                        'bcoop12','bcoop13','bcoop23','bcoop123',
#                                        'ucoop12','ucoop13','ucoop23','ucoop123']);

# 2DS np.array(['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1',
#            'f2', 'f12', 'kb1', 'ku1', 'kb2', 'ku2', 'bcoop', 'ucoop'])

# MDS np.array(['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm','fd', 'Km', 'Kd', 'kbm', 'kbd']);


for jj in range(4):

    for k in range(len(lists[jj])):
        mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
        myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
        myparams = np.loadtxt(file + namefiletosavedata + '_parms_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')

        # TODO replace numbers

        Kdim = myparams[20]/myparams[24]+myparams[21]/myparams[25]+myparams[22]/myparams[26]+myparams[23]/myparams[27]

        temp = mymeans_wt[vartoplot]
        if temp!=0:
            print(np.log10(Kdim))
            cm = plt.cm.get_cmap('rainbow')
            ll = ax.scatter(mymeans_wt[vartoplot], myvars_wt[vartoplot] / mymeans_wt[vartoplot], c=np.log10(Kdim), cmap=cm,vmin=np.log10(.05/20),vmax=np.log10(20/.05),  marker="o",s=5)

        ax.set_xlim([1,maxnumberofmolecules])
        ax.set_ylim([0,1000])

        ax.set_title(mymodel+'DS LrpB')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('dimer conc.')

cbar=fig.colorbar(ll, cax=axleg, label='cooperativity')
cbar.ax.set_yticklabels([str(10**(-2)),str(10**(-1)),'0','10',str(10**2)])
plt.tight_layout()

fig.savefig('noise_LrpB_coop.pdf')

