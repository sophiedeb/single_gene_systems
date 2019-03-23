"""plot_bifurcation.py: Make bifurcation plot."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from matplotlib import gridspec
from SGS import *
import matplotlib as mpl

mpl.rc('font',family='Calibri')

def cm2inch(value):
    return value/2.54

# ----------------------------------------------------------------------------------------------
# INPUT PARAMETERS
# ----------------------------------------------------------------------------------------------

sgs = DDS
system = System.RANDOM
concept = Concept.OSCILLATION

index = 2

# Path to save figure.
folder = foldername(sgs, system, concept)
figure_name = folder + "bifurcations/%d.pdf" % index

# ----------------------------------------------------------------------------------------------
# RECOVER PARAMETERS
# ----------------------------------------------------------------------------------------------

data = pd.read_csv(folder + "variables.csv", index_col=0)

params = data.loc[index].to_dict()

# Set fontsize.
mpl.rcParams['legend.fontsize'] = 12

# ----------------------------------------------------------------------------------------------
# RECOVER BIFURCATION DATA
# ----------------------------------------------------------------------------------------------

bif_data = pd.read_csv(folder + "bifurcations/bif-complete2.csv", index_col=0)

# -------------------------------------------------------------------------------
# GENERATE GRID FOR BIFURCATION DATA AND RECOVER LABELS
# -------------------------------------------------------------------------------

# Make grid for bifurcation plots.
if sgs == MDS:
    nrows = 4
    ncols = 5
    figsize = (8, 6)
    hr = [5, 5, 5, 1]
elif sgs == DDS:
    nrows = 5
    ncols = 4
    figsize = (cm2inch(18), cm2inch(18))
    hr = [5, 5, 5, 5, 1]
elif sgs == DDDS:
    nrows = 5
    ncols = 5
    figsize = (9, 7)
    hr = [5, 5, 5, 5, 1]

figev = plt.figure(figsize=figsize, tight_layout=True)

gs = gridspec.GridSpec(nrows,ncols, height_ratios=hr, hspace = 0.5, wspace=0.05)
axesev = [figev.add_subplot(gs[0])]
for i in range(1,(nrows-1)*ncols):
    axesev += [figev.add_subplot(gs[i], sharey=axesev[0])]

# Recover labels.
var, varlabel = parameter_list_latex(sgs)


# -------------------------------------------------------------------------------
# lowest row of grid for legends (period values, physiological range indication)
# -------------------------------------------------------------------------------

axesev += [figev.add_subplot(gs[(nrows-1)*ncols:(nrows-1)*ncols+int(ncols/2)])]
axesev += [figev.add_subplot(gs[(nrows-1)*ncols+int(ncols/2):])]

if len(varlabel)+1<nrows*ncols:
    for i in range(len(varlabel), len(axesev)-2):
        axesev[i].axis('off')

# -------------------------------------------------------------------------------
# color map of period
# -------------------------------------------------------------------------------

cmap = mpl.cm.get_cmap('inferno')

# -------------------------------------------------------------------------------
# ranges for amplitude and period
# -------------------------------------------------------------------------------

ampmax = np.nanmax(bif_data[bif_data.columns[pd.Series(bif_data.columns).str.endswith('max')]])
permin = np.nanmin(bif_data[bif_data.columns[pd.Series(bif_data.columns).str.endswith('per')]])
permax = np.nanmax(bif_data[bif_data.columns[pd.Series(bif_data.columns).str.endswith('per')]])

# -------------------------------------------------------------------------------
# fill grid with plots
# -------------------------------------------------------------------------------

for i, ax in enumerate(axesev):
    if i < len(varlabel):
        physmin, physmax = physiological_range(var[i])
        ax.set_title(varlabel[i], fontsize=9, y=0.75, x=0.05, ha='left')
        ax.set_xticks([])
        ax.axvline(params[var[i]], color='blue')
        ax.axvspan(physmin, physmax, alpha=0.2, color='red')
        #ax.set_xlim([0.2 * physmin, 2 * physmax])
        ax.set_xscale('log')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='both', which='major', labelsize=8)
    else:
        ax.axis('off')


for i, ax in enumerate(axesev):

    # -------------------------------------------------------------------------------
    # last cell of grid is filled with legend of period
    # -------------------------------------------------------------------------------

    if ax == axesev[-1]:
        N = 300
        gradient = np.linspace(0, 1, N)
        gradient = np.vstack((gradient, gradient))
        ax.imshow(gradient, cmap=cmap, aspect=N/30)
        #for j in range(len(scale)):
            #ax.plot(scale[j], 0, color=cmap((scale[j] - mins) / (maxs - mins)), marker='.')
        ax.set_yticks([])
        ax.set_xticks(np.array([0,0.25,0.5,0.75,1.0])*N)
        ax.margins(0.05)
        #ax.set_xticklabels(np.array([0, 0.25, 0.5, 0.75, 1.0])*(permax-permin)+permin)
        labels = [int(i) for i in np.array([0, 0.25, 0.5, 0.75, 1.0])*(permax-permin)+permin]
        ax.set_xticklabels(labels)
        ax.set_xlabel(r'Period (min)') #, fontsize=18)
        #ax.tick_params(labelsize=18)

    # -------------------------------------------------------------------------------
    # last but one cell is legend with physiological range indication
    # -------------------------------------------------------------------------------

    elif ax == axesev[-2]:
        h1, l1 = axesev[-3].get_legend_handles_labels()
        ax.axis('off')
        ax.legend(h1, l1, ncol=4, loc=10, frameon=False) #, fontsize=18)

    # -------------------------------------------------------------------------------
    # first cells, bifurcation diagrams of all different parameters
    # -------------------------------------------------------------------------------

    elif i < len(varlabel):
        # add title with parametername
        ax.set_title(varlabel[i], y=0.7, x=0.05, ha='left') # fontsize=18)

        # add shaded area denoting physiological ranges
        physmin = physiological_range(var[i])[0]
        physmax = physiological_range(var[i])[1]
        ax.axvspan(physmin, physmax, alpha=0.2, color='grey', label='physiological range')

        # get bifurcation data for the parameter
        x = 10**bif_data.loc[var[i], bif_data.columns[pd.Series(bif_data.columns).str.endswith('log10(par)')]].values
        bifmin = bif_data.loc[var[i], bif_data.columns[pd.Series(bif_data.columns).str.endswith('min')]].values
        bifmax = bif_data.loc[var[i], bif_data.columns[pd.Series(bif_data.columns).str.endswith('max')]].values
        period = bif_data.loc[var[i], bif_data.columns[pd.Series(bif_data.columns).str.endswith('per')]].values

        # plot the bifurcation diagram
        ax.scatter(x, bifmin, color=cmap((period-permin)/(permax-permin)), marker='.')
        ax.scatter(x, bifmax, color=cmap((period-permin)/(permax-permin)), marker='.')

        # denote the value of the parameter around which bifurcation diagram is made by black dotted line
        ax.axvline(x=params[var[i]], linestyle='dotted', color='k')

        # plotting specifics
        ax.set_ylim([-ampmax*0.1, ampmax])
        ax.set_xscale('log')
        ax.set_xlim([np.nanmin(x)-(x[1]-x[0]),np.nanmax(x)+(x[1]-x[0])])

        #ax.set_xticks([])
        ax.axis('on')

        # set labels
        if i == 0:
            ax.set_ylabel('Amplitude', rotation=0) #, fontsize=18)
            ax.yaxis.set_label_coords(-0.02, 1.05)
        if i%ncols!=0:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.tick_params(axis='both', which='major') #, labelsize=18)

# -------------------------------------------------------------------------------
# save figure
# -------------------------------------------------------------------------------

figev.savefig(figure_name)
plt.show()