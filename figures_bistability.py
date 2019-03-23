"""figures_bistability.py: Functions to analyze data and plot histograms, solutions in parameter space
and count for table."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from SGS import *
import pandas as pd
from scipy.spatial import ConvexHull
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib import rc

#from matplotlib.ticker import FormatStrFormatter
#plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%g'))

colorp = ['b', 'r', 'g', 'orange']

def min_max_values_3DS_physiological_range(v1, v2):
    ''' Returns the boundaries of the physiological ranges for 3DS of the A-F parameters.'''

    F = np.logspace(np.log10(3e-4), np.log10(3e4), 10)
    if v1== 'F' and v2== 'C':
        return F, 1e-3*F, F, 1e2*F
    elif v2 == 'F' and v1 == 'C':
        return 1e-3 * F, F, 1e2 * F, F
    else:
        Emin = 0.05/20*np.array([1e-8 + 2 * 1e-4 * (f - 2e-4) if f - 2e-4 <= 1e4
                else (1e-4 * (f - 1e4 - 1e-4) + 1e4 * (f - 1e4 - 1e-4) + 1) if f - 1e-4 - 1e4 <= 1e4
                else ((f - 2e4) ** 2 + 1e4 * (f - 2e4) + 1e8) for f in F])
        Emax = np.array([f ** 2 / 3 * 20/0.05 for f in F])
        E = np.logspace(np.log10(min(Emin)), np.log10(max(Emax)), 10)
        if v1 == 'F' and v2 == 'E':
            return F, Emin, F, Emax
        elif v2 == 'F' and v1 == 'E':
            return Emin, F, Emax, F
        elif v1 == 'C' and v2 == 'E':
            return np.append(1e-3*F[0], 1e2*F), np.append(Emin[0], Emin), np.append(1e-3*F,1e2*F[-1]), np.append(Emax, Emax[-1])
        elif v2 == 'C' and v1 == 'E':
            return np.append(Emin[0], Emin), np.append(1e-3*F[0], 1e2*F), np.append(Emax, Emax[-1]), np.append(1e-3*F,1e2*F[-1])
        elif v1 == 'C' and v2 == 'B':
            return np.append(1e-3*F[0], np.append(1e-3*F[0], 1e2*F)), np.append(1e-3*Emin[0], 1e-3*np.append(Emin[0], Emin)), \
                   np.append(np.append(1e-3*F,1e2*F[-1]), 1e2*F[-1]), np.append(1e2*np.append(Emax, Emax[-1]), 1e2*Emax[-1])
        elif v2 == 'C' and v1 == 'B':
            return np.append(1e-3 * Emin[0],1e-3 * np.append(Emin[0],Emin)), np.append(1e-3 * F[0], np.append(1e-3 * F[0], 1e2 * F)), \
                   np.append(1e2 * np.append(Emax, Emax[-1]), 1e2 * Emax[-1]), np.append(np.append(1e-3 * F, 1e2 * F[-1]), 1e2 * F[-1])
        elif v1== 'E' and v2== 'B':
            return E, 1e-3*E, E, 1e2*E
        elif v2 == 'E' and v1 == 'B':
            return 1e-3 * E, E, 1e2 * E, E
        elif v1== 'F' and v2== 'B':
            return F, 1e-3*Emin, F, 1e2*Emax
        elif v2 == 'F' and v1 == 'B':
            return 1e-3 * Emin, F, 1e2 * Emax, F
        else:
            Dmin = np.zeros(0)
            Dmax = np.zeros(0)
            for f in F:
                _, minD, _, _, _, _, _, _, maxD, _, _, _, _, _, _ = np.loadtxt('bistability/3DS/hull/%.3E-fig.txt' % f, unpack=True)
                Dmin = np.append(Dmin, min(minD))
                Dmax = np.append(Dmax, max(maxD))
            if v1 == 'F' and v2 == 'D':
                return F, Dmin, F, Dmax
            elif v2 == 'F' and v1 == 'D':
                return Dmin, F, Dmax, F
            elif v1 == 'F' and v2 == 'A':
                return F, 1e-3*Dmin, F, 1e2*Dmax
            elif v2 == 'F' and v1 == 'A':
                return 1e-3*Dmin, F, 1e2*Dmax, F
            elif v1 == 'C' and v2 == 'A':
                return np.append(1e-3*F[0], 1e2*F), 1e-3*np.append(Dmin[0], Dmin), \
                       np.append(1e-3*F, 1e2*F[-1]), 1e2*np.append(Dmax, Dmax[-1])
            elif v2 == 'C' and v1 == 'A':
                return 1e-3*np.append(Dmin[0], Dmin), np.append(1e-3*F[0], 1e2*F), \
                       1e2 * np.append(Dmax, Dmax[-1]), np.append(1e-3*F, 1e2*F[-1])
            elif v1 == 'C' and v2 == 'D':
                return np.append(1e-3 * F[0], 1e2 * F), np.append(Dmin[0], Dmin), \
                       np.append(1e-3 * F, 1e2 * F[-1]), np.append(Dmax, Dmax[-1])
            elif v2 == 'C' and v1 == 'D':
                return np.append(Dmin[0], Dmin), np.append(1e-3 * F[0], 1e2 * F), \
                       np.append(Dmax, Dmax[-1]), np.append(1e-3 * F, 1e2 * F[-1])

            Emin = np.zeros(0)
            Emax = np.zeros(0)
            Dmin = np.zeros(0)
            Dmax = np.zeros(0)
            for f in F:
                es, minDs, _, _, _, _, _,_, maxDs, _, _, _, _, _, _ = np.loadtxt('bistability/3DS/hull/%.3E-fig.txt' % f, unpack=True)
                for e, minD, maxD in zip(es, minDs, maxDs):
                    if len(Emin) < 2 or min(Emin) > e  or e > max(Emin):
                        Emin = np.append(Emin, e)
                        Emax = np.append(Emax, e)
                        Dmin = np.append(Dmin, minD)
                        Dmax = np.append(Dmax, maxD)
                    else:
                        i = int(np.where(Emin>e)[0][0])
                        yy = 10**(np.log10(Dmin[i-1]) + (np.log10(Dmin[i]) - np.log10(Dmin[i-1]))*(np.log10(e) - np.log10(Emin[i-1]))/\
                                                   (np.log10(Emin[i]) - np.log10(Emin[i-1])))
                        if minD < yy:
                            Emin = np.append(Emin, e)
                            Dmin = np.append(Dmin, minD)

                        i = int(np.where(Emax < e)[0][0])
                        yy = 10**(np.log10(Dmax[i - 1]) + (np.log10(Dmax[i]) - np.log10(Dmax[i - 1])) * (
                        np.log10(e) - np.log10(Emax[i - 1])) / (np.log10(Emax[i]) - np.log10(Emax[i - 1])))
                        if maxD > yy:
                            Emax = np.append(Emax,e)
                            Dmax = np.append(Dmax, maxD)
                s = sorted(zip(Emin, Dmin))
                Emin, Dmin = map(list, zip(*s))
                s = sorted(zip(Emax, Dmax))
                Emax, Dmax = map(list, zip(*s))
            todel = []
            for i in range(1,len(Emax)-1):
                if Dmax[i] < Dmax[i-1]:
                    todel += [i]
            Emax = np.delete(Emax, todel)
            Dmax = np.delete(Dmax, todel)

            todel = []
            for i in range(1, len(Emin) - 1):
                if Dmin[i] > Dmin[i+1]:
                    todel += [i]
            Emin = np.delete(Emin, todel)
            Dmin = np.delete(Dmin, todel)

            if v1 == 'E' and v2 == 'D':
                return Emin, Dmin, Emax, Dmax
            elif v2 == 'E' and v1 == 'D':
                return Dmin, Emin, Dmax, Emax
            elif v1 == 'A' and v2 == 'E':
                return np.append(1e-3 * Dmax[0], 1e2 * Dmax), np.append(Emax[0], Emax), \
                       np.append(1e-3 * Dmin, 1e2 * Dmin[-1]), np.append(Emin, Emin[-1])
            elif v2 == 'A' and v1 == 'E':
                return np.append(Emax[0], Emax), np.append(1e-3 * Dmax[0], 1e2 * Dmax),\
                       np.append(Emin, Emin[-1]), np.append(1e-3 * Dmin, 1e2 * Dmin[-1])
            elif v1 == 'B' and v2 == 'D':
                return np.append(1e-3*Emin[0], 1e2*Emin), np.append(Dmin[0], Dmin), \
                       np.append(1e-3*Emax,1e2*Emax[-1]), np.append(Dmax, Dmax[-1])
            elif v2 == 'B' and v1 == 'D':
                return np.append(Dmin[0], Dmin), np.append(1e-3*Emin[0], 1e2*Emin), \
                       np.append(Dmax, Dmax[-1]), np.append(1e-3 * Emax, 1e2 * Emax[-1])
            elif v1 == 'B' and v2 == 'A':
                return np.append(1e-3*Emin[0], 1e2*Emin), 1e-3*np.append(Dmin[0], Dmin), \
                       np.append(1e-3*Emax,1e2*Emax[-1]), 1e2*np.append(Dmax, Dmax[-1])
            elif v2 == 'B' and v1 == 'A':
                return 1e-3 * np.append(Dmin[0], Dmin), np.append(1e-3 * Emin[0], 1e2 * Emin), \
                       1e2 * np.append(Dmax,Dmax[-1]), np.append(1e-3 * Emax, 1e2 * Emax[-1]),
            elif v1 == 'D' and v2 == 'A':
                D = np.logspace(np.log10(min(Dmin)), np.log10(max(Dmax)), 10)
                return D, 1e-3*D, D, 1e2*D
            elif v2 == 'D' and v1 == 'A':
                D = np.logspace(np.log10(min(Dmin)), np.log10(max(Dmax)), 10)
                return 1e-3*D, D, 1e2*D, D
    return np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)

def histograms(ax, os, f, zoom=False, xlabel=True):
    Nbins = 10

    data = pd.read_csv(f, na_values='NAN')
    data = data[np.logical_and(data['nonmono'] != 4, np.log10(data['T'])>0)]
    data = data[~np.isnan(np.log10(data['T']))]

    T = data['T']

    if os == DDS:
        #data = pd.read_csv(f)
        #_, T, _ = np.loadtxt(ft, unpack=True)
        #_, nm = np.loadtxt(fnm, unpack=True)
        x1 = 2.7; x2 = 4.2
    elif os == MDS:
        #_, _, _, _, _, T = np.loadtxt(ft, unpack=True)
        #_, nm = np.loadtxt(fnm, unpack=True)
        x1 = 3.0
        x2 = 4.25
    elif os == DDDS:
        #data = np.loadtxt(ft, unpack=True)
        #if len(data[0])>3:
        #    data = data[:3, :]
        #_, T, _ = data
        #_, nm = np.loadtxt(fnm, unpack=True)
        x1 = 2.2
        x2 = 4.1

    if not zoom:
        T = np.log10(T)
    else:
        data = data[T<600]
        T = T[T<600]

    _, bins = np.histogram(T, bins=Nbins)
    ax.hist(T[data['nonmono'] == 0], bins=bins, histtype=u'step', edgecolor=colorp[0], label="monotonic")
    ax.hist(T[data['nonmono'] == 1], bins=bins, histtype=u'step', edgecolor=colorp[1], label="non-monotonic 1")
    ax.hist(T[data['nonmono'] == 2], bins=bins, histtype=u'step', edgecolor=colorp[2], label="non-monotonic 2")
    ax.hist(T[data['nonmono'] == 3], bins=bins, histtype=u'step', edgecolor=colorp[3], label="non-monotonic 3")

    if not zoom:
        ax.set_xlim([x1, x2])
    ax.set_yscale('log')
    if xlabel and not zoom:
        ax.set_xlabel(r'$\Delta t^* [\log_{10}$(min)]')
    elif xlabel and zoom:
        ax.set_xlabel(r'$\Delta t^* [min]$')

def induction_time_histograms_figure():
    ''' Plot the 4 histograms with the distribution of induction time.'''

    # Set font
    font = {'family': 'Arial', 'size': 10}
    rc('font', **font)

    # Make figure and grid.
    fig = plt.figure(figsize=(5.2,2.7), tight_layout=True)
    gs = gridspec.GridSpec(2,4, wspace=0.2, hspace=0.7, height_ratios=[3,1])
    gs.tight_layout(fig)

    wl = -0.085 #width distance label
    hl = 1.2 #height distance label

    # Add all subfigures (MDS, 2DS, 3DS, 3DS zoomed)
    ax0 = fig.add_subplot(gs[2])
    ax0.text(wl, hl, 'C. 3DS', transform=ax0.transAxes, va='top', ha='left')
    histograms(ax0, DDDS, 'bistability/3DS/H500-600.csv',
               xlabel=False)
    ax0.tick_params(labelleft=False)

    ax = fig.add_subplot(gs[0],sharey=ax0)
    ax.text(wl, hl, 'A. MDS', transform=ax.transAxes, va='top', ha='left')
    histograms(ax, MDS, 'bistability/MDS/H500-600.csv', xlabel = False)

    ax = fig.add_subplot(gs[1], sharey=ax)
    ax.text(wl, hl, 'B. 2DS', transform=ax.transAxes, va='top', ha='left')
    histograms(ax, DDS, 'bistability/2DS/H500-600.csv')
    ax.tick_params(labelleft=False)

    ax = fig.add_subplot(gs[3], sharey=ax0)
    ax.text(wl, hl, 'D. 3DS, zoomed', transform=ax.transAxes, va='top', ha='left')
    histograms(ax, DDDS, 'bistability/3DS/H500-600.csv', zoom=True)
    ax.tick_params(labelleft=False)

    # Plot legend.
    axleg = fig.add_subplot(gs[1,:])
    legend_elements = [Line2D([0], [0], color=colorp[0], lw=2, label='monotonic'),
                        Line2D([0], [0], color=colorp[1], lw=2, label='non-monotonic 1'),
                        Line2D([0], [0], color=colorp[2], lw=2, label='non-monotonic 2'),
                        Line2D([0], [0], color=colorp[3], lw=2, label='non-monotonic 3'), ]

    axleg.legend(handles=legend_elements, loc='upper center', fontsize=10, handlelength=1, ncol=2, frameon=False)
    axleg.axis('off')

    # Save and show figure
    #fig.savefig('Figures/bistabparhists-zoom.eps')

    plt.show()

induction_time_histograms_figure()

def bistability_in_parameter_space_figure():

    # Ticks in integer format instead of float
    def my_formatter(x, pos):
        val_str = '{:g}'.format(x)
        if val_str == '0.0':
            val_str = '0'
        return val_str

    from matplotlib.ticker import FuncFormatter
    major_formatter = FuncFormatter(my_formatter)

    # Set font sizes
    font = {'family': 'Arial', 'size': 10}

    rc('font', **font)
    rc('xtick', labelsize=8)
    rc('ytick', labelsize=8)

    # Define grids for figures (gs upper grid for MDS and 2DS, gs2 for lower grid for 3DS)
    wspace = 0.05
    hspace = 0.05
    gs = gridspec.GridSpec(13, 11, height_ratios=[1,1,1,1,1,1.5,1,1,1,1,1,1,1], width_ratios=[1,1,1,1,1,1.3,1,1,1,1,1], wspace = wspace, hspace = hspace)
    gs2 = gridspec.GridSpec(13, 9, height_ratios=[1,1,1,1,1,1.5,1,1,1,1,1,1,1], width_ratios=[1,1,1,1,1,1,1,0.3,3], wspace = wspace, hspace = hspace)

    # Label parameters (font size and distances)
    width_distance_label_subfigure = -0.05
    height_distance_label_subfigure = 1.5
    coordinate_y_labels = -0.6
    label_size = 8

    # Make figure
    fig = plt.figure(figsize=(6, 8))

    # Function that returns set of axes
    def make_axes(labels, gs, ox=0, oy=0):
        ax = [[0 for _  in range(len(labels))] for _ in range(len(labels))]
        diag_ax = [0 for _  in range(len(labels))]

        for i in range(len(labels)):
            for j in range(len(labels)):
                if i == 0 and j == 0:
                    ax[i][j] = fig.add_subplot(gs[i + oy, j+ox])
                    diag_ax[i] = ax[i][j]._make_twin_axes(sharex=ax[i][j], frameon=False)
                    diag_ax[i].axis('off')
                elif j == 0:
                    ax[i][j] = fig.add_subplot(gs[i + oy, j+ox], sharex=ax[0][0])
                elif i == 0:
                    ax[i][j] = fig.add_subplot(gs[i + oy, j+ox], sharey=ax[0][0])
                elif i == j:
                    ax[i][j] = fig.add_subplot(gs[i + oy, j+ox], sharey=diag_ax[0], sharex=ax[0][j])
                    diag_ax[i] = ax[i][j]._make_twin_axes(sharex=ax[i][j], sharey=diag_ax[0], frameon=False)
                    diag_ax[i].axis('off')
                else:
                    ax[i][j] = fig.add_subplot(gs[i + oy, j+ox], sharex=ax[0][j], sharey=ax[i][0])

                ax[i][j].tick_params(axis='both', which='major', length=3)
                if i != len(labels) - 1:
                    plt.setp(ax[i][j].get_xticklabels(), visible=False)
                    ax[i][j].tick_params(axis='x', which='major', labelsize=label_size, length=3)
                if j != 0:
                    plt.setp(ax[i][j].get_yticklabels(), visible=False)
                    ax[i][j].tick_params(axis='y', which='major', labelsize=label_size, length=3)

        for i in range(len(labels)):
            ax[-1][i].set_xlabel(labels[i])
            ax[i][0].set_ylabel(labels[i])
            ax[i][0].get_yaxis().set_label_coords(coordinate_y_labels, 0.5)
            ax[i][0].yaxis.set_major_formatter(major_formatter)
            ax[-1][i].xaxis.set_major_formatter(major_formatter)

        return ax, diag_ax

    # Function that fills the axes with the data
    def fill_axes_data(data, ax, diag_ax, parname, parameter_ranges):
        for i in range(len(parname)):
            for j in range(len(parname)):
                if i == j:
                    vals = []
                    for k in range(4):
                        # Attempt to get data for this level, allowing for empty
                        try:
                            vals.append(np.log10(data[parname[i]][data['nonmono'] == k]))
                        except KeyError:
                            vals.append(np.array([]))

                    diag_ax[i].hist(vals, color=colorp, histtype="barstacked", bins=7, range=parameter_ranges[i])
                elif i>j:
                    ax[i][j].scatter(np.log10(data[parname[j]][data['nonmono'] == 1]),
                                     np.log10(data[parname[i]][data['nonmono'] == 1]), c=colorp[1], s=1)
                    ax[i][j].scatter(np.log10(data[parname[j]][data['nonmono'] == 2]),
                                     np.log10(data[parname[i]][data['nonmono'] == 2]), c=colorp[2], s=1)
                    ax[i][j].scatter(np.log10(data[parname[j]][data['nonmono'] == 3]),
                                     np.log10(data[parname[i]][data['nonmono'] == 3]), c=colorp[3], s=1)
                else:
                    ax[i][j].scatter(np.log10(data[parname[j]][data['nonmono'] == 0]),
                                     np.log10(data[parname[i]][data['nonmono'] == 0]), c=colorp[0], s=1)

    # Plot different subfigures

    plotMDS = True
    plot2DS = True
    plot3DS = True

    if plotMDS:
        # Read data
        f = 'bistability/MDS/H500-600.csv'

        data = pd.read_csv(f, na_values='NAN')

        # Define parameters, and ranges
        parname = ['Kd', 'Km', 'fd', 'fm', 'T']
        labels = [r'$K_d^*$', r'$K_m^*$', '$f_d^*$', '$f_m^*$', r'$\Delta t^*$']

        parameter_ranges = [(-7,3), (-10,3), (-3,2), (-3,2), (3,4)]

        # Make the axes, write the label and plot the data
        ax, diag_ax = make_axes(labels, gs)

        ax[0][0].text(width_distance_label_subfigure, height_distance_label_subfigure, 'A. MDS',
                      transform=ax[0][0].transAxes, va='top', ha='right')

        fill_axes_data(data, ax, diag_ax, parname, parameter_ranges)

        # Contour the physiological parameter space
        for i in range(len(labels)):
            for j in range(i+1,len(labels)):
                if not(i == len(labels) -1 or j == len(labels) - 1):
                    ax[i][j].plot([parameter_ranges[j][0], parameter_ranges[j][0], parameter_ranges[j][1], parameter_ranges[j][1], parameter_ranges[j][0]], [parameter_ranges[i][1], parameter_ranges[i][0], parameter_ranges[i][0], parameter_ranges[i][1], parameter_ranges[i][1]], 'k:', linewidth=0.5)
                    ax[j][i].plot([parameter_ranges[i][0], parameter_ranges[i][0], parameter_ranges[i][1], parameter_ranges[i][1], parameter_ranges[i][0]], [parameter_ranges[j][1], parameter_ranges[j][0], parameter_ranges[j][0], parameter_ranges[j][1], parameter_ranges[j][1]], 'k:', linewidth=0.5)

    if plot2DS:
        # Read data
        f = 'bistability/2DS/H500-600.csv'

        data = pd.read_csv(f, na_values='NAN')

        # Define parameters, and ranges
        parname = ['A', 'B', 'C', 'D', 'T']
        labels = [r'$A^*$', r'$B^*$', '$C^*$', '$D^*$', r'$\Delta t^*$']
        parameter_ranges = [(-15, 7), (-8, 8), (-12, 5), (-5, 5), (2.7, 4)]

        # Make the axes, write the label and plot the data
        ax, diag_ax = make_axes(labels, gs, ox=6)
        ax[0][0].text(width_distance_label_subfigure, height_distance_label_subfigure, 'B. 2DS',
                      transform=ax[0][0].transAxes, va='top', ha='right')

        fill_axes_data(data, ax, diag_ax, parname, parameter_ranges)

        # Contour the physiological parameter space
        for i in range(len(labels)):
            for j in range(i+1,len(labels)):
                if not(i == len(labels) -1 or j == len(labels) - 1):
                    p = np.hstack((np.array([np.log10(data[parname[j]])]).T, np.array([np.log10(data[parname[i]])]).T))
                    hull = ConvexHull(p)
                    for simplex in hull.simplices:
                        print(simplex)
                        ax[i][j].plot(p[simplex, 0], p[simplex, 1], 'k:', linewidth=0.5)
                        ax[j][i].plot(p[simplex, 1], p[simplex, 0], 'k:', linewidth=0.5)

    if plot3DS:
        # Read data
        f = 'bistability/3DS/H500-600.csv'

        data = pd.read_csv(f, na_values='NAN')
        data = data[np.logical_and(np.logical_and(data['nonmono'] < 4, data['T'] >= -1), ~np.isnan(data['T']))]

        # Define parameters, and ranges
        parname = ['A', 'B', 'C', 'D', 'E', 'F', 'T']
        labels = [r'$A^*$', r'$B^*$', '$C^*$', '$D^*$', '$E^*$', '$F^*$', r'$\Delta t^*$']
        parameter_ranges = [(-29.9, 29.9), (-14.9, 14.9), (-8, 8), (-24.9, 24.9), (-14.9, 14.9), (-4.5, 4.5), (2.55, 3.9)]

        # Make the axes, write the label and plot the data
        ax, diag_ax = make_axes(labels, gs2, oy=6)
        ax[0][0].text(width_distance_label_subfigure, height_distance_label_subfigure, 'C. 3DS',
                      transform=ax[0][0].transAxes, va='top', ha='right')
        fill_axes_data(data, ax, diag_ax, parname, parameter_ranges)

        # Contour the physiological parameter space
        for i in range(len(labels)):
            for j in range(i+1, len(labels)):
                if not(i == len(labels) -1 or j == len(labels) - 1):
                    x1, y1, x2, y2 = min_max_values_3DS_physiological_range(parname[j], parname[i])
                    ax[i][j].plot(np.log10(x1), np.log10(y1), 'k:', linewidth=0.5)
                    ax[i][j].plot(np.log10(x2), np.log10(y2), 'k:', linewidth=0.5)
                    if len(x1) > 0:
                        ax[i][j].plot([np.log10(x1[0]), np.log10(x2[0])], [np.log10(y1[0]), np.log10(y2[0])], 'k:',
                                linewidth=0.5)
                        ax[i][j].plot([np.log10(x1[-1]), np.log10(x2[-1])], [np.log10(y1[-1]), np.log10(y2[-1])], 'k:',
                                linewidth=0.5)

                    ax[j][i].plot(np.log10(y1), np.log10(x1), 'k:', linewidth=0.5)
                    ax[j][i].plot(np.log10(y2), np.log10(x2), 'k:', linewidth=0.5)
                    if len(x1) > 0:
                        ax[j][i].plot([np.log10(y1[0]), np.log10(y2[0])], [np.log10(x1[0]), np.log10(x2[0])], 'k:',
                                      linewidth=0.5)
                        ax[j][i].plot([np.log10(y1[-1]), np.log10(y2[-1])], [np.log10(x1[-1]), np.log10(x2[-1])],
                                      'k:', linewidth=0.5)

    # Plot legend.
    axleg = fig.add_subplot(gs2[6:,7:], frameon=False)

    legend_elements = [
        Line2D([0], [0], color='w', lw=2, label='monotonic', marker='o', markerfacecolor=colorp[0], markersize=3),
        Line2D([0], [0], color='w', lw=2, label='non-monotonic 1', marker='o', markerfacecolor=colorp[1], markersize=3),
        Line2D([0], [0], color='w', lw=2, label='non-monotonic 2', marker='o', markerfacecolor=colorp[2], markersize=3),
        Line2D([0], [0], color='w', lw=2, label='non-monotonic 3', marker='o', markerfacecolor=colorp[3], markersize=3)]

    axleg.legend(handles=legend_elements, loc='center', fontsize=10, handlelength=0.5, ncol=1, frameon=False)
    axleg.axis('off')

    # Save and show figure
    #gs.tight_layout(fig)

    #fig.savefig('Figures/bistabpars.png')
    plt.show()

#bistability_in_parameter_space_figure()

def count_bistability_solutions():
    systems = ['MDS', '2DS', '3DS']

    for system in systems:
        f = 'bistability/%s/H500-600.csv' % system
        data = pd.read_csv(f, na_values='NAN')
        data = data[np.logical_and(data['nonmono'] != 4, np.log10(data['T']) > 0)]
        data.dropna(inplace=True)

        print(system, data.groupby('nonmono').count()['R'])
        print("total %s" % system, len(data))

#count_bistability_solutions()

def time_versus_score_figure(f):
    data = pd.read_csv(f)

    plt.figure()
    plt.plot(data['R'], data['T'], '.')
    plt.xlabel('Score R')
    plt.ylabel('Time [min]')
    plt.xscale('log')
    plt.yscale('log')

#time_versus_score_figure('bistability/MDS/H500-600.csv')
#plt.show()

