"""figures_oscillation.py: Functions to analyze data and plot histograms, solutions in parameter space."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from matplotlib import gridspec
from SGS import *
from pylab import *
from matplotlib.lines import Line2D
import os.path
import matplotlib as mpl
import pandas as pd

colorp = ['b', 'r', 'g', 'orange']


class Project(Enum):
    PAPER = 1
    SLIDE = 2

def cm2inch(value):
    return value / 2.54


def logarithmic_parameter_volume(data, sgs, system=System.RANDOM):
    """ Returns logarithmic volume."""

    # Determine all variables that participate in the volume.
    if sgs == DDDS:
        if system == System.SSLRPB:
            var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12',
                   'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'bcoop12', 'bcoop23', 'bcoop13', 'bcoop123']
        elif system == System.RANDOM:
            var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12',
                   'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop23', 'bcoop13',
                   'bcoop123', 'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']
    elif sgs == DDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12', 'kb1',
               'kb2', 'ku1', 'ku2', 'bcoop', 'ucoop']
    elif sgs == MDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm',
               'kbd']

    # Determine logarithmic volume.
    edges = data['log10_right_bound'].loc[var] - data['log10_left_bound'].loc[var]
    vol = np.prod(edges)

    if vol < 0:
        print("ERROR: VOLUME SMALLER THAN 0")
    return vol


def ax_histogram_variable(ax, data, var, varlabel):
    """" Plot histogram of variable. """

    # Set physiological range.
    physmin = physiological_range(var)[0]
    physmax = physiological_range(var)[1]
    ax.axvspan(physmin, physmax, alpha=0.2, color='red')
    ax.set_xlim([0.2 * physmin, 2 * physmax])

    # Set title, ticks, axes labels.
    ax.set_title(varlabel, fontsize=9, y=0.75, x=0.05, ha='left')
    ax.set_xticks([])
    ax.set_xscale('log')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='both', which='major', labelsize=8)

    # Get data.
    if isinstance(var, str):
        d = data[var]
    else:
        d = np.zeros(0)
        for j in var:
            d = np.append(d, data[:, j])

    # Plot histogram.
    bins = np.logspace(np.log10(physmin), np.log10(physmax), 20)
    ax.hist(d, bins=bins)


def histograms_variables_figure(data, os, system=System.RANDOM):
    """ Make a figure with the histograms of all variables of the oscillators."""

    # Set the list of variables.
    if os == DDDS:
        if system == System.SSLRPB:
            var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12',
                   'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop23', 'bcoop13',
                   'bcoop123', 'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']
            varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                        r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_1$', r'$f_2$', r'$f_3$',
                        r'$f_{12}$', r'$f_{13}$', r'$f_{23}$', r'$f_{123}$', r'$k_{b1}$', r'$k_{b2}$', r'$k_{b3}$',
                        r'$k_{u1}$', r'$k_{u2}$', r'$k_{u3}$', r'co$_{b12}$', r'co$_{b23}$', r'co$_{b13}$',
                        r'co$_{b123}$', r'co$_{u12}$', r'co$_{u23}$', r'co$_{u13}$', r'co$_{u123}$']

            nrows = 6
            ncols = 5
            figsize = (9, 4)
        elif system == System.RANDOM:
            var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', ['f1', 'f2', 'f3'],
                   ['f12', 'f13', 'f23'], 'f123', ['kb1', 'kb2', 'kb3'], ['ku1', 'ku2', 'ku3'],
                   ['bcoop12', 'bcoop23', 'bcoop13'], 'bcoop123', ['ucoop12', 'ucoop23', 'ucoop13'], 'ucoop123']
            varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                        r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_i$', r'$f_{ij}$',
                        r'$f_{123}$', r'$k_{bi}$', r'$k_{ui}$', r'co$_{bij}$', r'co$_{b123}$', r'co$_{uij}$',
                        r'co$_{u123}$']

            nrows = 4
            ncols = 4
            figsize = (9, 4)
    elif os == DDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', ['f1', 'f2'], 'f12',
               ['kb1', 'kb2'], ['ku1', 'ku2'], 'bcoop', 'ucoop']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_i$', r'$f_{12}$', r'$k_{b}$',
                    r'$k_{u}$', r'co$_{b}$', r'co$_{u}$']

        nrows = 4
        ncols = 4
        figsize = (7, 7)
    elif os == MDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm',
               'kbd']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_m$', r'$f_d$', r'$K_{m}$',
                    r'$K_{d}$', r'$k_{bm}$', r'$k_{bd}$']

        nrows = 4
        ncols = 4
        figsize = (7, 7)

    # Set the figure and subfigures.
    fig = plt.figure('parameters', figsize=figsize, tight_layout=True)
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.3, wspace=0.3)
    axesev = [fig.add_subplot(i) for i in gs]

    # Draw the histograms in the different subfigures.
    for i, ax in enumerate(axesev):
        if i < len(varlabel):
            ax_histogram_variable(ax, data, var[i], varlabel[i])
        else:
            ax.axis('off')


def plot_bifurcations(fname, data, folder, os, system=System.RANDOM):
    """ Overview of ranges of selection of solutions."""

    if os == DDDS:
        if system == System.SSLRPB:
            var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12',
                   'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop23', 'bcoop13',
                   'bcoop123', 'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']
            varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                        r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_1$', r'$f_2$', r'$f_3$',
                        r'$f_{12}$', r'$f_{13}$', r'$f_{23}$', r'$f_{123}$', r'$k_{b1}$', r'$k_{b2}$', r'$k_{b3}$',
                        r'$k_{u1}$', r'$k_{u2}$', r'$k_{u3}$', r'co$_{b12}$', r'co$_{b23}$', r'co$_{b13}$',
                        r'co$_{b123}$', r'co$_{u12}$', r'co$_{u23}$', r'co$_{u13}$', r'co$_{u123}$']

            nrows = 5
            ncols = 6
            figsize = (cm2inch(15), cm2inch(15))
        elif system == System.RANDOM:
            var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1',
                   'f12', 'f123', 'kb1', 'ku1',
                   'bcoop12', 'bcoop123', 'ucoop12', 'ucoop123']
            varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                        r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_i$', r'$f_{ij}$',
                        r'$f_{123}$', r'$k_{bi}$', r'$k_{ui}$', r'co$_{bij}$', r'co$_{b123}$', r'co$_{uij}$',
                        r'co$_{u123}$']

            nrows = 3
            ncols = 6
            figsize = (9, 4)
            figsize = (cm2inch(15), cm2inch(10))

    elif os == DDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f12',
               'kb1', 'ku1', 'bcoop', 'ucoop']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_i$', r'$f_{12}$', r'$k_{b}$',
                    r'$k_{u}$', r'co$_{b}$', r'co$_{u}$']

        nrows = 4
        ncols = 4
        figsize = (cm2inch(15), cm2inch(15))

    elif os == MDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm',
               'kbd']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_m$', r'$f_d$', r'$K_{m}$',
                    r'$K_{d}$', r'$k_{bm}$', r'$k_{bd}$']

        nrows = 4
        ncols = 4
        figsize = (cm2inch(15), cm2inch(15))
    N = 20 # number of shown solutions

    # Make figure and grid for subfigures.
    fig = plt.figure('bifurcations', figsize=figsize, tight_layout=True)
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.3, wspace=0.3)
    axesev = [fig.add_subplot(i) for i in gs]
    axesev += [fig.add_subplot(gs[-2:])]


    for i, ax in enumerate(axesev):
        if i < len(varlabel):
            # Draw physiological ranges on figure.
            physmin = physiological_range(var[i])[0]
            physmax = physiological_range(var[i])[1]
            ax.axvspan(physmin, physmax, alpha=0.2, color='red')

            # Title, ticks, labels.
            ax.set_title(varlabel[i], fontsize=9, y=0.75, x=0.05, ha='left')
            ax.set_xticks([])
            ax.set_xlim([0.2 * physmin, 2 * physmax])
            ax.set_xscale('log')
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='both', which='major', labelsize=8)
            ax.set_yticks([])



        elif i == len(axesev) - 1:
            # Legend.
            legend_elements = [Line2D([0], [0], color='b', lw=2, label='monotonic'),
                               Line2D([0], [0], color='r', lw=2, label='non-monotonic 1'),
                               Line2D([0], [0], color='g', lw=2, label='non-monotonic 2'),
                               Line2D([0], [0], color='y', lw=2, label='non-monotonic 3'), ]

            ax.legend(handles=legend_elements, loc='center', fontsize=8, handlelength=1)
            ax.axis('off')
        else:
            ax.axis('off')

    # Draw ranges for different solutions.
    jj = 0
    colors = ['blue', 'red', 'g', 'y']
    for i in range(len(colors)):
        df = data[data['monotonicity_type'] == i]
        indices = df[:min(len(df), int(N / len(data) * len(df)))].index.tolist()

        for idx in indices:
            data_idx = pd.read_csv(folder + "%d.csv" % idx, index_col=0)
            for j in range(len(var)):
                left_bound = 10 ** data_idx['log10_left_bound'].loc[var[j]]
                right_bound = 10 ** data_idx['log10_right_bound'].loc[var[j]]

                axesev[j].plot([left_bound, right_bound], [jj, jj], c=colors[i])
                jj += 1

    if(fname != 0):
        fig.savefig(fname)
    plt.show()


def histograms_volume_distributions_figure(fname, data, os, folder, title, system=System.RANDOM):
    """ Figure of the distributions of the logarithmic volumes separated by monotonicity type."""

    # Set figure, subfigure grid and title.
    fig = plt.figure(figsize=(3, 4.5), tight_layout=True)
    gs = gridspec.GridSpec(5, 2, hspace=0., height_ratios=[2, 5, 5, 5, 5], width_ratios=[0, 5])
    fig.suptitle(title, x=0.1, y=0.95, fontsize=11)

    # Calculate the volumes for every monotonicity type.
    vols = [[], [], [], []]
    for i in range(4):
        vol = np.zeros(0)
        for j in data[data['monotonicity_type'] == i].index.tolist():
            vol = np.append(vol,
                            logarithmic_parameter_volume(pd.read_csv(folder + "%d.csv" % j, index_col=0), os, system))
        vols[i] = vol

    # Plot the histograms with labels.
    a = np.log10(np.nanmin([np.nanmin(i) if len(i) > 0 else np.nan for i in vols]))
    b = np.log10(np.nanmax([np.nanmax(i) if len(i) > 0 else np.nan for i in vols]))
    bins = np.logspace(a - 0.1 * (b - a), b + 0.1 * (b - a), 15)
    nonmonolabel = ['Monotonic', 'Non-monotonic 1', 'Non-monotonic 2', 'Non-monotonic 3']
    for i in range(4):
        if i == 0:
            ax0 = fig.add_subplot(gs[i + 1, 1])
            ax = ax0
        else:
            ax = fig.add_subplot(gs[i + 1, 1], sharex=ax0)
        ax.set_title(nonmonolabel[i], fontsize=11, y=0.65, x=0.05, ha='left')

        # Plot the histogram and print the results.
        if len(vols[i]) > 0:
            ax.hist(vols[i], bins=bins)
            print(nonmonolabel[i], "median, mean, std:", '%.3E' % np.median(vols[i]), '%.3E' % np.mean(vols[i]),
                  '%.3E' % np.std(vols[i]))
            if not os == MDS:
                ax.set_yscale('log')
            ax.set_xscale('log')
        else:
            ax.set_yticks([])
            print(nonmonolabel[i], "no results")

        if i < 3:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel('Logarithmic parameter volume')
    vols = np.concatenate(vols).ravel()  # flatten vols list
    print('total', '%.3E' % np.median(vols), '%.3E' % np.mean(vols), '%.3E' % np.std(vols))

    # Set ylabel.
    fig.text(0.02, 0.5, 'Number of solutions', va='center', rotation='vertical')

    # gs.tight_layout(fig, h_pad=0) #, rect=[0.5, 0, 1, 1], h_pad=0.5)

    # Save and plot figure.
    fig.savefig(fname)
    plt.show()


def histograms_combined_parameters_figure(data, os):
    """ Make figure of the parameters that are combinations of other parameters (Ki, Kdim, omega)."""

    # Determine the figure structure and variables given the oscillator.
    if os == DDDS:
        nrows = 4
        ncols = 2
        figsize = (9, 9)

        var = [('kb1', 'ku1'), ('kb2', 'ku2'), ('kb3', 'ku3'), ('alphaass', 'alphadiss'), ('bcoop12', 'ucoop12'),
               ('bcoop13', 'ucoop13'), ('bcoop23', 'ucoop23'), ('bcoop123', 'ucoop123')]
        varlabel = ['$K_1$', '$K_2$', '$K_3$', '$K_{dim}$', r'$\omega_{12}$', r'$\omega_{13}$', r'$\omega_{23}$',
                    r'$\omega_{123}$']
    elif os == DDS:
        nrows = 2
        ncols = 2
        figsize = (4, 4)

        var = [('kb1', 'ku1'), ('kb2', 'ku2'), ('alphaass', 'alphadiss'), ('bcoop', 'ucoop')]
        varlabel = ['$K_1$', '$K_2$', '$K_{dim}$', r'$\omega$']
    elif os == MDS:
        nrows = 2
        ncols = 2
        figsize = (4, 4)

        var = [('kbm', 'kum'), ('kbd', 'kud'), ('alphaass', 'alphadiss')]
        varlabel = ['$k_{um}$', '$k_{ud}$', '$K_{dim}$']

    # Make the figure and grid for subfigures.
    fig = plt.figure('parameters 2', figsize=figsize, tight_layout=True)
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.3, wspace=0.3)

    # Plot histograms in all subfigures.
    for i in range(len(var)):
        ax = fig.add_subplot(gs[i])

        # Color physiological range.
        physmin = physiological_range(var[i][0])[0] / physiological_range(var[i][1])[1]
        physmax = physiological_range(var[i][0])[1] / physiological_range(var[i][1])[0]
        ax.axvspan(physmin, physmax, alpha=0.2, color='red')

        # Set title, ticks and labels.
        ax.set_title(varlabel[i], fontsize=9, y=0.75, x=0.05, ha='left')
        ax.set_xticks([])
        ax.set_xscale('log')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='both', which='major', labelsize=8)

        # Plot histogram.
        a = data[var[i][0]] / data[var[i][1]]
        bins = np.logspace(np.log10(min(a)), np.log10(max(a)), 20)
        ax.hist(a, bins=bins)

    plt.show()


def histograms_amplitudes_figure(data, os):
    """ Plot the histograms of the amplitudes of all variables of the oscillator."""

    # Set the variables depending on the oscillator and size of figure.
    if os == DDDS:
        var = ['DNA0', 'DNA1', 'DNA2', 'DNA3', 'DNA12', 'DNA13', 'DNA23', 'DNA123', 'mRNA', 'm', 'd']
        varlabel = [r'$DNA_0$', r'$DNA_1$', r'$DNA_2$', r'$DNA_{3}$', r'$DNA_{12}$', r'$DNA_{13}$', r'$DNA_{23}$',
                    r'$DNA_{123}$', r'$mRNA$', r'$m$', r'$d$']

        nrows = 4
        ncols = 3
        figsize = (9, 6)
    elif os == DDS:
        var = ['DNA0', 'DNA1', 'DNA2', 'DNA12', 'mRNA', 'm', 'd']
        varlabel = [r'$DNA_0$', r'$DNA_1$', r'$DNA_2$', r'$DNA_{12}$', r'$mRNA$', r'$m$', r'$d$']

        nrows = 4
        ncols = 2
        figsize = (9, 4)
    elif os == MDS:
        var = ['DNA0', 'DNAm', 'DNAd', 'mRNA', 'm', 'd']
        varlabel = [r'$DNA_0$', r'$DNA_m$', r'$DNA_d$', r'$mRNA$', r'$m$', r'$d$']

        nrows = 4
        ncols = 2
        figsize = (9, 4)

    # Make figure and grid for subfigures.
    fig = plt.figure('amplitudes', figsize=figsize, tight_layout=True)
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.3, wspace=0.3)
    axesev = [fig.add_subplot(i) for i in gs]

    # Plot histogram in every subfigure.
    for i, ax in enumerate(axesev):
        if i < len(varlabel):
            # Title, ticks, labels.
            ax.set_title(varlabel[i], fontsize=9, y=0.75, x=0.05, ha='left')
            ax.set_xticks([])
            ax.set_yscale('log')
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='both', which='major', labelsize=8)

            # Plot histogram.
            amplitudes = data['%s_max' % var[i]] - data['%s_min' % var[i]]
            ax.hist(amplitudes)
        else:
            ax.axis('off')

    plt.show()


def histograms_mean_figure(data, os):
    """ Plot histogram of means of oscillation of all variables of the oscillator."""

    # Set the variable lists.
    if os == DDDS:
        nrows = 5
        ncols = 3
        varlabel = [r'$DNA_0$', r'$DNA_1$', r'$DNA_2$', r'$DNA_{3}$', r'$DNA_{12}$', r'$DNA_{13}$', r'$DNA_{23}$',
                    r'$DNA_{123}$', r'$mRNA$', r'$m$', r'$d$']
        var = ['DNA0', 'DNA1', 'DNA2', 'DNA3', 'DNA12', 'DNA13', 'DNA23', 'DNA123', 'mRNA', 'm', 'd']
        figsize = (15, 10)
    elif os == DDS:
        varlabel = [r'$DNA_0$', r'$DNA_1$', r'$DNA_2$', r'$DNA_{12}$', r'$mRNA$', r'$m$', r'$d$']
        var = ['DNA0', 'DNA1', 'DNA2', 'DNA12', 'mRNA', 'm', 'd']
        nrows = 4
        ncols = 3
        figsize = (7, 7)
    elif os == MDS:
        nrows = 4
        ncols = 3
        figsize = (7, 7)
        varlabel = [r'$DNA_0$', r'$DNA_m$', r'$DNA_d$', r'$mRNA$', r'$m$', r'$d$']
        var = ['DNA0', 'DNAm', 'DNAd', 'mRNA', 'm', 'd']
    # Variables of which we will compare the means.
    relative_data = [('m', 'mRNA'), ('d', 'm'), ('d', 'mRNA')]

    # Set figure and subfigures.
    fig = plt.figure('means', figsize=figsize, tight_layout=True)
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.8, wspace=0.3)
    axesev = [fig.add_subplot(i) for i in gs]

    # Plot histograms in the subplots.
    for i, ax in enumerate(axesev):
        if i < len(varlabel):
            # Title, labels, ticks
            ax.set_title(varlabel[i], fontsize=9, y=0.75, x=0.05, ha='left')
            ax.set_yscale('log')
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='both', which='major', labelsize=8)

            # Plot the histogram of the means.
            means = (data['%s_max' % var[i]] + data['%s_min' % var[i]]) / 2
            ax.hist(means)
        elif i < len(varlabel) + len(relative_data):
            # Title, labels, ticks.
            ax.set_title(relative_data[i - len(varlabel)][0] + '/' + relative_data[i - len(varlabel)][1], fontsize=9,
                         y=0.75, x=0.05, ha='left')
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(axis='both', which='major', labelsize=8)

            # Plot the histogram of the relative means.
            means = data[relative_data[i - len(varlabel)][0] + '_max'] / data[
                relative_data[i - len(varlabel)][1] + '_max']
            ax.hist(means, bins=np.logspace(np.log10(min(means)), np.log10(max(means)), 20))
        else:
            ax.axis('off')

    plt.show()


def protein_ranges_3d_figure(data):
    """ Make figure with the protein (mRNA, m, d) ranges of the oscillations in a 3D plot."""

    # Make figure.
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111, projection='3d')

    # Get data.
    minm = data['m_min'].tolist()
    maxm = data['m_max'].tolist()
    mind = data['d_min'].tolist()
    maxd = data['d_max'].tolist()
    minmRNA = data['mRNA_min'].tolist()
    maxmRNA = data['mRNA_max'].tolist()

    # Plot the ranges of all solutions.
    from mpl_toolkits import mplot3d
    c = ['blue' if nm == 0 else 'red' if nm == 1 else 'green' if nm == 2 else 'y' for nm in
         data['monotonicity_type'].tolist()]
    for i in range(len(maxm)):
        ax.plot3D([minm[i], maxm[i]], [minmRNA[i], maxmRNA[i]], [mind[i], maxd[i]], color=c[i])

    # Labels.
    ax.set_xlabel('m')
    ax.set_ylabel('mRNA')
    ax.set_zlabel('d')
    ax.invert_yaxis()

    # Show figure.
    plt.show()


def protein_ranges_2d_figure(fname, data):
    """ Make figure with the protein (m, d) ranges of the oscillations in monomer-dimer plane."""

    # Set figure.
    fig = plt.figure(figsize=(3, 3), tight_layout=True)
    ax = fig.add_subplot(111)

    # Get data.
    minm = data['m_min'].tolist()
    maxm = data['m_max'].tolist()
    mind = data['d_min'].tolist()
    maxd = data['d_max'].tolist()

    # Plot ranges.
    from mpl_toolkits import mplot3d
    c = ['blue' if nm == 0 else 'red' if nm == 1 else 'green' if nm == 2 else 'y' for nm in
         data['monotonicity_type'].tolist()]
    # ax.scatter(means[:,-2], means[:,-1], means[:,-3], c=c)
    for i in range(len(maxm)):
        ax.plot([minm[i], maxm[i]], [mind[i], maxd[i]], color=c[i])

    # Plot diagonal (m = d).
    ax.plot([0, min(max(maxm), max(maxd))], [0, min(max(maxm), max(maxd))], 'k:')

    # Labels.
    ax.set_xlabel('m')
    ax.set_ylabel('d')

    # Save and show figure.
    fig.savefig(fname)
    plt.show()


def histogram_period_figure(data):
    """ Make figure with histogram of the periods of the oscillators."""

    period = data['period']
    plt.figure('period')
    plt.hist(period, bins=np.linspace(min(period), max(period), 20))
    plt.yscale('log')

    plt.show()


def ax_amplitudes(ax, data):
    """ Draw the amplitudes in the monomer-dimer plane."""

    # Get data.
    maxm = data['m_max'].tolist()
    maxd = data['d_max'].tolist()
    monotonicity_type = data['monotonicity_type'].tolist()

    # Draw dots.
    c = [colorp[0] if nm == 0 else colorp[1] if nm == 1 else colorp[2] if nm == 2 else colorp[3] for nm in
         monotonicity_type]
    ax.scatter(maxm, maxd, color=c, s=3)

    # Plot the diagonal.
    ax.plot([0, max(max(maxm), max(maxd))], [0, max(max(maxm), max(maxd))], 'k:')

    # Set labels.
    ax.set_xlabel('$\max(m)$')
    ax.set_ylabel('$\max(d)$')

    # Set axes equal (square).
    ax.axis('equal')

    # Set limits.
    # ax.set_xlim(-.05*max(max(maxm),max(maxd)), 1.05*max(max(maxm),max(maxd)))
    ax.set_ylim(-.05 * max(max(maxm), max(maxd)), 0.6 * max(max(maxm), max(maxd)))


def amplitudes_figure(project=Project.PAPER):
    """ Create figure with all amplitudes."""

    # Set font and font size.
    font = {'family': 'Arial', 'size': 10}
    mpl.rc('font', **font)
    mpl.rc('legend', handlelength=1)

    # Make figure and grid for subfigures.
    fig = plt.figure(figsize=(5.2, cm2inch(8)), tight_layout=True)
    gs = gridspec.GridSpec(2, 3, height_ratios=[5, 5], width_ratios=[2, 2, 1], hspace=1.5, wspace=0.8)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    axleg = fig.add_subplot(gs[:, 2])

    os = ['MDS', '2DS', '3DS', 'SsLrpB']

    # Define titles of subfigures.
    if project == Project.PAPER:
        t = ['A. MDS', 'B. 2DS', 'C. 3DS', 'D. Ss-LrpB']
        wl = -0.125  # width distance label
        hl = 1.2  # height distance label
    elif project == Project.SLIDE:
        t = ['MDS', '2DS', '3DS', 'SsLrpB']
        # wl = 0.9  # width distance label
        # hl = 0.85  # height distance label

    # Generate all subfigures.
    for i, ax in enumerate(axes):
        if i < len(os):
            # Load data.
            data = pd.read_csv('oscillation/' + os[i] + '/variables.csv')
            data = data[~np.isnan(data['period'])]

            # Draw amplitudes.
            ax_amplitudes(ax, data)

            # Set title.
            if project == Project.PAPER:
                ax.text(wl, hl, t[i], transform=ax.transAxes, va='top', ha='left')
            elif project.SLIDE:
                ax.text(wl, hl, t[i], transform=ax.transAxes, ha='right', size=14)
        else:
            ax.axis('off')

    # Draw legend.
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=colorp[0], lw=2, label='monotonic'),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=colorp[1], lw=2,
                              label='non-monotonic 1'),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=colorp[2], lw=2,
                              label='non-monotonic 2'),
                       Line2D([0], [0], marker='o', color='w', markerfacecolor=colorp[3], lw=2,
                              label='non-monotonic 3'), ]
    axleg.legend(handles=legend_elements, loc='center', handlelength=1, ncol=1, frameon=False)
    axleg.axis('off')

    # Align the labels.
    fig.align_labels(axes)

    # Save and show the figure.
    gs.tight_layout(fig)
    fig.savefig('Figures/amplitudes4.eps')
    plt.show()


def ax_parameter_ranges(ax, data, folder, os, system=System.RANDOM, concept=Concept.OSCILLATION):
    """ Plot the mean parameter and mean parameter range of the set of oscillators in the axis."""

    if os == DDDS:
        if system == System.SSLRPB:
            var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12',
                   'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop23', 'bcoop13',
                   'bcoop123', 'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']
            varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                        r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_1$', r'$f_2$', r'$f_3$',
                        r'$f_{12}$', r'$f_{13}$', r'$f_{23}$', r'$f_{123}$', r'$k_{b1}$', r'$k_{b2}$', r'$k_{b3}$',
                        r'$k_{u1}$', r'$k_{u2}$', r'$k_{u3}$', r'co$_{b12}$', r'co$_{b23}$', r'co$_{b13}$',
                        r'co$_{b123}$', r'co$_{u12}$', r'co$_{u23}$', r'co$_{u13}$', r'co$_{u123}$']
        elif system == System.RANDOM:
            var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', ['f1', 'f2', 'f3'],
                   ['f12', 'f13', 'f23'], 'f123', ['kb1', 'kb2', 'kb3'], ['ku1', 'ku2', 'ku3'],
                   ['bcoop12', 'bcoop23', 'bcoop13'], 'bcoop123', ['ucoop12', 'ucoop23', 'ucoop13'], 'ucoop123']
            varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                        r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_i$', r'$f_{ij}$',
                        r'$f_{123}$', r'$k_{bi}$', r'$k_{ui}$', r'co$_{bij}$', r'co$_{b123}$', r'co$_{uij}$',
                        r'co$_{u123}$']
    elif os == DDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', ['f1', 'f2'], 'f12',
               ['kb1', 'kb2'], ['ku1', 'ku2'], 'bcoop', 'ucoop']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_i$', r'$f_{12}$', r'$k_{b}$',
                    r'$k_{u}$', r'co$_{b}$', r'co$_{u}$']
    elif os == MDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm',
               'kbd']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_m$', r'$f_d$', r'$K_{m}$',
                    r'$K_{d}$', r'$k_{bm}$', r'$k_{bd}$']

    # Plot the physiological range between 0 and 1.
    ax.axvspan(0, 1, color='gainsboro', zorder=1)

    # Plot range of every parameter.
    for i in range(len(varlabel)):
        print(var[i])
        # Get the bounds and widths of all oscillators.
        bounds = np.zeros(0)
        widths = np.zeros(0)
        if isinstance(var[i], str):
            # Get physiological ranges of parameter.
            physmin = physiological_range(var[i])[0]
            physmax = physiological_range(var[i])[1]

            if concept == Concept.OSCILLATION:
                data = data[~np.isnan(data['period'])]

            # Set indices
            indices = data.index.tolist()

            for j in indices:
                data_j = pd.read_csv(folder + "%d.csv" % j, index_col=0)
                bounds = np.append(bounds, data_j['log10_left_bound'].loc[var[i]])
                bounds = np.append(bounds, data_j['log10_right_bound'].loc[var[i]])
                widths = np.append(widths,
                                   data_j['log10_right_bound'].loc[var[i]] - data_j['log10_left_bound'].loc[var[i]])
        else:
            # Get physiological ranges of parameter.
            physmin = physiological_range(var[i][0])[0]
            physmax = physiological_range(var[i][0])[1]

            if concept == Concept.OSCILLATION:
                indices = data[~np.isnan(data['period'])].index.tolist()

            for j in indices:
                data_j = pd.read_csv(folder + "%d.csv" % j, index_col=0)
                for v in var:
                    bounds = np.append(bounds, data_j['log10_left_bound'].loc[v])
                    bounds = np.append(bounds, data_j['log10_right_bound'].loc[v])
                    widths = np.append(widths, data_j['log10_right_bound'].loc[v] - data_j['log10_left_bound'].loc[v])

        # Calculate the mean and mean width, rescale between 0 and 1 for the physiological range.
        m = (np.mean(np.array(bounds)) - np.log10(physmin)) / (np.log10(physmax) - np.log10(physmin))
        w = (np.mean(np.array(widths))) / (np.log10(physmax) - np.log10(physmin))

        # Plot the mean and the range (alternatingly in black and grey for clearness).
        if i % 2 == 0:
            c = 'k'
        else:
            c = 'gray'
        ax.scatter([m], [len(varlabel) - i], s=5, c=c, zorder=2)
        ax.plot([m - w / 2, m + w / 2], [len(varlabel) - i, len(varlabel) - i], c=c, zorder=3)

    # Set limits, ticks.
    ax.set_xlim([-0.2, 1.2])
    ax.set_xticks([])
    ax.set_yticks(range(1, len(varlabel) + 1))
    ax.set_yticklabels(reversed(varlabel))
    for i, ytick in zip(range(len(varlabel)), ax.get_yticklabels()):
        if (len(varlabel) - i) % 2 == 0:
            ytick.set_color('gray')
        else:
            ytick.set_color('k')
        ytick.set_fontsize(10)


def parameter_ranges_figure(concept = Concept.OSCILLATION):
    """ Make figure of the parameter ranges for all systems (MDS, 2DS, 3DS, SsLrpB)."""

    # Set font and font size.
    font = {'family': 'Arial', 'size': 10}
    mpl.rc('font', **font)
    mpl.rc('legend', handlelength=1)

    # Plot the ranges for every system in a separate subfigure.
    if concept == Concept.OSCILLATION:
        oss = [(MDS, System.RANDOM, 'MDS'), (DDS, System.RANDOM, '2DS'), (DDDS, System.RANDOM, '3DS'),
               (DDDS, System.SSLRPB, 'SsLrpB')]
        t = ['A. MDS', 'B. 2DS', 'C. 3DS', 'D. Ss-LrpB']
    elif concept == Concept.BISTABILITY:
        oss = [(MDS, System.RANDOM, 'MDS'), (DDS, System.RANDOM, '2DS'), (DDDS, System.RANDOM, '3DS')]
        t = ['A. MDS', 'B. 2DS', 'C. 3DS']

    # Set figure and grid for subfigures.
    fig = plt.figure(figsize=(5.2, 6), tight_layout=True)
    gs = gridspec.GridSpec(1, len(oss), wspace=0.75)
    axesev = [fig.add_subplot(i) for i in gs]

    # Set folder name to oscillation/bistability
    if concept == Concept.OSCILLATION:
        f = 'oscillation/'
    elif concept == Concept.BISTABILITY:
        f = 'bistability/'

    for i, ax in enumerate(axesev):
        print(t[i])
        data = pd.read_csv(f + oss[i][2] + '/variables.csv', index_col=0)

        if concept == Concept.OSCILLATION:
            data = data[~np.isnan(data['period'])]

        folder = f + oss[i][2] + '/bifurcations/points-'
        ax_parameter_ranges(ax, data, folder, oss[i][0], oss[i][1], concept)
        ax.set_title(t[i], loc='left', fontsize=10)
    gs.tight_layout(fig)

    # Save and show figure.
    fig.savefig('Figures/parameterbifurcation-%s.eps' % f[:-1])
    plt.show()


def main():
    # add_oscillation_data("oscillation/MDS/variables.csv", MDS)
    # add_selection_monotonicity_data("oscillation/SsLrpB/variables.csv", DDDS)
    # count_solutions("oscillation/SsLrpB/variables.csv")

    #file = 'SsLrpB'
    #os = DDDS
    # fname = 'oscillation/' + file + '/volumedistribution.pdf'
    #folder = 'oscillation/' + file + '/bifurcations/'
    #title = file
    #system = System.SSLRPB if file == 'SsLrpB' else System.RANDOM
    #data = pd.read_csv("oscillation/" + file + "/variables.csv", index_col=0)
    #data = data[~np.isnan(data['period'])]

    # volume_distributions_figure(fname, data, os, folder, title, system)

    #file = '3DS'
    #os = DDDS
    #fname = 0 # 'oscillation/' + file + '/parameterbifurcation.pdf'
    #folder = 'oscillation/' + file + '/bifurcations/points-'
    #system = System.SSLRPB if file == 'SsLrpB' else System.RANDOM
    #data = pd.read_csv("oscillation/" + file + "/variables.csv", index_col=0)
    #data = data[~np.isnan(data['period'])]

    # plot_bifurcations(fname, data, folder, os, system)

    # combined_parameters_figure(data, os)

    # histograms_amplitudes_figure(data, os)

    # histograms_mean_figure(data, os)

    # protein_ranges_3d_figure(data)

    # fname = 'oscillation/' + file + '/protein.pdf'
    # protein_ranges_2d_figure(fname, data)

    # histogram_period_figure(data)

    # amplitudes_figure()

    # Fig 4 of paper
    #parameter_ranges_figure(concept=Concept.BISTABILITY)

    return


if __name__ == "__main__":
    main()
