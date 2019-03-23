"""plot_SGS.py: Visualization of SGS, timeseries, steady state, eigenvalues, parameters in parameter space,
solutions in phase space with nullclines."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

import matplotlib as mpl
from SGS import *
from matplotlib import gridspec
import pandas as pd

# ----------------------------------------------------------------------------------------------
# INPUT PARAMETERS
# ----------------------------------------------------------------------------------------------

sgs = DDDS

stochsim = False  # choose whether or not to perform stochastic simulation

filename = "oscillation/SsLrpB/variables.csv"
index = 1250

save = 0  # choose whether you want to save the figures
data = pd.read_csv(filename, index_col=0)

mpl.rcParams['legend.fontsize'] = 12

DNAcopynumber = 1.0

# ----------------------------------------------------------------------------------------------
# RECOVER PARAMETERS AND MAKE SGS OBJECT
# ----------------------------------------------------------------------------------------------

# Recover parameters.
params = data.loc[index].to_dict()

print("l", np.log10(params['f123']), np.log10(params['f13']))
# Some parameters are scaled with the number of cells (DNA copy number).
scaledparams = ['kb1', 'kb2', 'kb3', 'kbd', 'kbm', 'alphaass']

for i in scaledparams:
    if i in params:
        params[i] /= DNAcopynumber

ddds = sgs(params)

ddds.DNAtot = DNAcopynumber

# Print parameters and if they are compatible with the SsLrpB system.
print("%s with parameters:" % sgs.__name__, ddds.parameter_line())
if sgs == DDDS:
    print("DDDS is SsLrpB compatible", ddds.SsLrpB_compatible())

# ----------------------------------------------------------------------------------------------
# MAKE STEADY-STATE FIGURE
# ----------------------------------------------------------------------------------------------

ssfig = plt.figure('Steady-state ratios and transcription rate', figsize=(8, 5), tight_layout=True)
ss = ddds.steady_state()
print("Steady state", ss)
ax1m = ssfig.add_subplot(2, 2, 1)
ddds.plot_quasi_steady_state(ax1m)
ax2m = ssfig.add_subplot(2, 2, 3)
ddds.plot_rate_balance(ax2m)

ssx = [ss[a] for a in sgs.ALL_VARS]

ax1d = ssfig.add_subplot(2, 2, 2)
ddds.plot_quasi_steady_state(ax1d, polymer=Polymer.DIMER)
ax2d = ssfig.add_subplot(2, 2, 4)
ddds.plot_rate_balance(ax2d, polymer=Polymer.DIMER, maxx=1000)

# ----------------------------------------------------------------------------------------------
# MAKE NULLCLINES FIGURE
# ----------------------------------------------------------------------------------------------

print('Nullclines')
ncfig = plt.figure('Nullclines', figsize=(4, 4), tight_layout=True)
ax3 = ncfig.add_subplot(2, 1, 1)
ddds.plot_nullclines(ax3)
ax4 = ncfig.add_subplot(2, 1, 2)
ddds.plot_nullclines(ax4, polymer=Polymer.DIMER)

# ----------------------------------------------------------------------------------------------
# MAKE DETERMINISTIC TIMELAPSE FIGURE
# ----------------------------------------------------------------------------------------------

t_f = 5000
dt = 0.5

hf = {'m': lambda t: 0.0, 'd': lambda t: 0.0, 'mRNA': lambda t: 1.00, 'DNA0': lambda t: DNAcopynumber,
      'DNA1': lambda t: 0.0, 'DNA2': lambda t: 0.0, 'DNA12': lambda t: 0.0, 'DNA3': lambda t: 0.0,
      'DNA13': lambda t: 0.0, 'DNA23': lambda t: 0.0, 'DNA123': lambda t: 0.0, }

print('timelapse')
if len(ss['d']) == 1:
    fig = plt.figure('timelapse prot', figsize=(4, 5), tight_layout=True)
    gs = mpl.gridspec.GridSpec(3, 1, height_ratios=[4, 4, 1])  # 1 -> space for legend
    ax5 = fig.add_subplot(gs[0])
    ax6 = fig.add_subplot(gs[1])

    hf = {'m': lambda t: 0, 'd': lambda t: 0, 'mRNA': lambda t: 1.0, 'DNA0': lambda t: 1.0, 'DNA1': lambda t: 0,
          'DNA2': lambda t: 0, 'DNA3': lambda t: 0, 'DNA12': lambda t: 0, 'DNA13': lambda t: 0, 'DNA23': lambda t: 0,
          'DNA123': lambda t: 0, 'DNAm' : lambda t: 0, 'DNAd' : lambda t: 0}

    dde = ddds.time_series(t_f, dt, hf=hf)
    ddds.plot_timelapse(ax5, dde, t_f, dt, DNA=False, axncm = ax3, axncd = ax4)
    ddds.plot_timelapse(ax6, dde, t_f, dt, DNA=True)
    if save:
        fig.savefig('timelapse-3DO.pdf')
        # np.savetxt('oscillation/SsLrpB/stochastic/367/deterministicSsLrpB-1603-1.txt', dde.sol['d'][10000:])
else:
    fig = plt.figure('timelapse prot', figsize=(5, 3.5), tight_layout=True)

    starts = [0, 15, 30, 45, 50, 58, 100, 150, 200, 250, 300, 350, 400, 500, 600]
    for i in starts:
        hf = {'m': lambda t: np.sqrt(i * (ddds.alphadiss + ddds.gammad) / ddds.alphaass), 'd': lambda t: i,
              'mRNA': lambda t: max(0, int(ss['mRNA'][1]) - 1), 'DNA0': lambda t: 1, 'DNA1': lambda t: 0,
              'DNA2': lambda t: 0, 'DNA3': lambda t: 0, 'DNA12': lambda t: 0, 'DNA13': lambda t: 0,
              'DNA23': lambda t: 0, 'DNA123': lambda t: 0, 'DNAm' : lambda t: 0, 'DNAd' : lambda t: 0}

        dde = ddds.time_series(t_f, dt, hf=hf)
        ax = fig.add_subplot(111)
        ddds.plot_timelapse(ax, dde, t_f, dt, DNA=False, legend=(i == 0), only='d') #, axncm = ax3, axncd = ax4)
        # plt.gca().set_color_cycle(None)

        # fig.savefig("bistability/bistability2.pdf")

# ----------------------------------------------------------------------------------------------
# MAKE STOCHASTIC TIMELAPSE FIGURE
# ----------------------------------------------------------------------------------------------

t_f = 50
dt = 0.1
if stochsim:
    figtsstoch = plt.figure('timelapse stochastic', figsize=(4, 5), tight_layout=True)
    for i in [0]: #range(len(ss['m'])):
        hf = np.array([round(ss['DNA0'][i]), round(ss['DNA1'][i]), round(ss['DNA2'][i]), round(ss['DNA3'][i]),
                        round(ss['DNA12'][i]), round(ss['DNA23'][i]), round(ss['DNA13'][i]), round(ss['mRNA'][i]),
                        round(ss['m'][i]), 60])

        axa = figtsstoch.add_subplot(2, 1, 1)
        dde2 = ddds.stochastic_time_series(t_f, dt, hf=hf)
        ddds.plot_timelapse(axa, dde2, t_f, dt, DNA=False, axncm=ax3, axncd=ax4)
        axb = figtsstoch.add_subplot(2, 1, 2)
        ddds.plot_timelapse(axb, dde2, t_f, dt, DNA=True)
        np.set_printoptions(linewidth=3000, suppress=True)
        #np.savetxt('oscillation/SsLrpB/stochastic275e-20.txt', dde2.sol['d'])

        # ----------------------------------------------------------------------------------------------
        # MAKE DETERMINISTIC-STOCHASTIC COMBINATION TIMELAPSE FIGURE
        # ----------------------------------------------------------------------------------------------

        if False:
            fig = plt.figure(figsize=(5, 3), tight_layout=True)
            gs = gridspec.GridSpec(2, 1, hspace=0)
            ax = fig.add_subplot(gs[0])
            ddds.plot_timelapse(ax, dde, t_f, dt, t_i=500)
            ax = fig.add_subplot(gs[1], sharex=ax)
            ddds.plot_timelapse(ax, dde2, t_f, dt, legend=False)
            plt.savefig('oscillation/SsLrpB/timelapse306-1.pdf')

print('Eigenvalues')
figev = plt.figure('Eigenvalues', figsize=(3, 2), tight_layout=True)
ax7 = figev.add_subplot(1, 1, 1)
a = ddds.eigenvalues(ax7)

# ----------------------------------------------------------------------------------------------
# MAKE PARAMETER IN RANGE FIGURE
# ----------------------------------------------------------------------------------------------

if True:
    var, varlabel = parameter_list_latex(sgs)

    if sgs == MDS:
        nrows = 4
        ncols = 4
        figsize = (6, 6)
    elif sgs == DDS:
        nrows = 4
        ncols = 4
        figsize = (6, 6)
    elif sgs == DDDS:
        nrows = 5
        ncols = 6
        figsize = (9, 4)

    fig = plt.figure('parameters', figsize=figsize, tight_layout=True)
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.3, wspace=0.3)

    axesev = [fig.add_subplot(i) for i in gs]

    for i, ax in enumerate(axesev):
        if i < len(varlabel):
            physmin, physmax = physiological_range(var[i])
            ax.set_title(varlabel[i], fontsize=9, y=0.75, x=0.05, ha='left')
            ax.set_xticks([])
            ax.axvline(params[var[i]], color='blue')
            ax.axvspan(physmin, physmax, alpha=0.2, color='red')
            ax.set_xlim([0.2 * physmin, 2 * physmax])
            ax.set_xscale('log')
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='both', which='major', labelsize=8)
        else:
            ax.axis('off')

hf = {'m': lambda t: 0.0, 'd': lambda t: 0.0, 'mRNA': lambda t: 1.00, 'DNA0': lambda t: 1.0, 'DNA1': lambda t: 0.0,
      'DNA2': lambda t: 0.0, 'DNA12': lambda t: 0.0, 'DNA3': lambda t: 0.0, 'DNA13': lambda t: 0.0,
      'DNA23': lambda t: 0.0, 'DNA123': lambda t: 0.0, }

t_f = 1500
dt = 0.1
print("eigenvalues", a)
print("is oscillatory", ddds.is_oscillatory())
print("is bistable", ddds.is_bistable())

#print("oscillation parameters", ddds.oscillation_parameters(ddds.ALL_VARS))

plt.show()
