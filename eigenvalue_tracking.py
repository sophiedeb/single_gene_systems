"""eigenvalue_tracking.py: Plot eigenvalues for small variations of parameters of solution."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

import matplotlib as mpl
from matplotlib import gridspec
from SGS import *

# ----------------------------------------------------------------------------------------------
# INPUT PARAMETERS
# ----------------------------------------------------------------------------------------------

oscillator = DDDS;

filename = "oscillation/3DS/variables.csv"
index = 45

save = 0  # choose whether you want to save the figures
data = pd.read_csv(filename, index_col=0)

mpl.rcParams['legend.fontsize'] = 12

DNAcopynumber = 1.0

# ----------------------------------------------------------------------------------------------
# RECOVER PARAMETERS AND MAKE DDDS OBJECT
# ----------------------------------------------------------------------------------------------

params = data.loc[index].to_dict()


"""
PARAMETERS
"""

if (oscillator == DDS):
    varname = ['kb1', 'ku1', 'kb2', 'ku2', 'beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'bcoop',
               'ucoop', 'phi0', 'fa', 'fr']
    var = [params[i] for i in varname];
    varlabel = [r'$k_{b1}$', r'$k_{u1}$', r'$k_{b2}$', r'$k_{u2}$', r'$\beta$', r'$\gamma_\mathrm{m}$',
                r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$', r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$',
                r'coop$_b$', r'coop$_u$', r'$\phi_0$', r'$f_a$', r'$f_r$']
    nrows = 4
    ncols = 4
    figsize = (6, 6)
elif (oscillator == MDS):
    varname = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fa', 'fr', 'kbm', 'kum',
                   'kbd', 'kud'];
    varlabel = varname[:]
    var = [params[i] for i in varname];
    nrows = 4
    ncols = 4
    figsize = (6, 6)
elif (oscillator == DDDS):
    varname = ['kb1', 'ku1', 'kb2', 'ku2', 'kb3', 'ku3', 'beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass',
               'alphadiss', 'bcoop12', 'ucoop12', 'bcoop13', 'ucoop13', 'bcoop23', 'ucoop23', 'bcoop123', 'ucoop123',
               'phi0', 'f1', 'f2', 'f3', 'f12', 'f23', 'f13', 'f123']
    var = [params[i] for i in varname];
    varlabel = [r'$k_{b1}$', r'$k_{u1}$', r'$k_{b2}$', r'$k_{u2}$', r'$k_{b3}$', r'$k_{u3}$', r'$\beta$',
            r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$', r'$\alpha_\mathrm{ass}$',
            r'$\alpha_\mathrm{diss}$', r'coop$_{b12}$', r'coop$_{u12}$', r'coop$_{b13}$', r'coop$_{u13}$',
            r'coop$_{b23}$', r'coop$_{u23}$', r'coop$_{b123}$', r'coop$_{u123}$', r'$\phi_0$', r'$f_1$', r'$f_2$', r'$f_3$',
            r'$f_{12}$', r'$f_{23}$', r'$f_{13}$', r'$f_{123}$']
    nrows = 5
    ncols = 6
    figsize=(9, 7)

mins = -1.5;
maxs = 1.5;
scale = np.linspace(mins, maxs, 20);
cmap = mpl.cm.get_cmap('coolwarm')

#axes = np.ones(len(var));
minev = np.ones(len(var));
maxev = np.ones(len(var));

figev = plt.figure(figsize=figsize, tight_layout=True)
#fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(7, 5))
#axes = axes.flatten()
gs = gridspec.GridSpec(nrows,ncols, hspace = 0, wspace=0)

axesev = [figev.add_subplot(i) for i in gs]

for i in range(len(var)):
    print(varlabel[i])
    axesev[i].set_title(varlabel[i], fontsize=9, y=0.6, x=0.05, ha='left')
    for j in range(len(scale)):
        params[varname[i]] = var[i] * 10 ** scale[j]
        do = oscillator(params)
        do.eigenvalues(ax=axesev[i], color=cmap((scale[j] - mins) / (maxs - mins)))

    params[varname[i]] = var[i]
    do = oscillator(params)
    do.eigenvalues(ax=axesev[i], color='k')

    # eig = eig.flatten();
    # minev[i] = min(np.real(eig));
    # maxev[i] = max(np.real(eig));

if(len(var)+1<nrows*ncols):
    for i in range(len(var), nrows*ncols-1):
        axesev[i].axis('off')

for ax in axesev:
    if(ax == axesev[-1]):
        N = 300
        gradient = np.linspace(0, 1, N)
        gradient = np.vstack((gradient, gradient))
        ax.imshow(gradient, cmap=cmap, aspect=N/10)
        #for j in range(len(scale)):
            #ax.plot(scale[j], 0, color=cmap((scale[j] - mins) / (maxs - mins)), marker='.')
        ax.set_yticks([])
        ax.set_xticks(np.array([0,0.25,0.5,0.75,1.0])*N)
        ax.margins(0.05)
        ax.set_xticklabels([-1,-0.5,0,0.5,1])
        #ax.set_xticklabels(mins+np.array([0, 0.25, 0.5, 0.75, 1.0])*(maxs-mins))
        ax.set_xlabel(r'$s$', fontsize=8)
        ax.tick_params(labelsize=8)
    else:
        # ax.set_xlim([1.05 * min(minev), 1.05 * max(maxev)]);
        ax.set_xlim([-0.15, 0.15])
        #ax.set_ylim([-0.15, 0.15])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')

plt.show()