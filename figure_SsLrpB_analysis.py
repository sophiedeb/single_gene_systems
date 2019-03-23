"""figure_SsLrpB_analysis.py: Calculate number of SsLrpB compatible solutions, plot figure with compatible region,
bistable region and oscillatory solutions."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from singlegenesystem import *
import matplotlib as mpl
import os.path
from collections import OrderedDict
import pandas as pd


class Project(Enum):
    PAPER = 1
    SLIDE = 2


def percentage_SsLrpB_compatible_solutions(fname):
    """ Save 3DS solutions compatible with SsLrpB.

    Returns matrix, every row contains values for different f13-f123 combinations:
    f13, f123, min value of C for which solution is compatible, max value of C with compatibility
    and number of compatible solutions out of 100"""

    # SsLrpB measured values.
    ct = 1e-3 / 2.4

    Kd1 = 73.5 * ct
    Kd2 = 0.7 * ct
    Kd3 = 49.1 * ct
    omega12 = 4.1
    omega13 = 7.9
    omega23 = 2.1
    omega123 = 3.1

    d = omega12 * omega23 * omega13 * omega123 * Kd1 * Kd2 * Kd3
    e = omega13 * Kd1 * Kd3 + omega12 * Kd1 * Kd2 + omega23 * Kd2 * Kd3
    f = Kd1 + Kd2 + Kd3

    # Free parameters.
    f123s = np.logspace(-3, 2, 51)
    Cs = np.logspace(-3 + np.log10(Kd1 + Kd3), 2 + np.log10(Kd1 + Kd3), 100)
    f13s = np.logspace(-3, 2, 51)

    # Arrays to store found solutions.
    goodf13 = np.zeros(0)
    goodf123 = np.zeros(0)
    cmin = np.zeros(0)
    cmax = np.zeros(0)
    N = np.zeros(0)

    # Loop over free parameters.
    for f123 in f123s:
        a = f123 * omega12 * omega23 * omega13 * omega123 * Kd1 * Kd2 * Kd3
        for f13 in f13s:
            b = f13 * omega13 * Kd1 * Kd3

            cs = np.zeros(0)
            for i in range(len(Cs)):
                c = Cs[i]

                xxx = np.linspace(0, 1000, 10000)

                fy = lambda x: (a * x ** 3 + b * x ** 2 + c * x + 1) / (d * x ** 3 + e * x ** 2 + f * x + 1)

                y = fy(xxx)

                if max(y) > 2.0 and min(y) < 0.5 and np.argmax(y) < np.argmin(y) and min(y) == y[-1]:
                    if False:
                        plt.figure()
                        plt.plot(xxx, y)
                        plt.show()
                    cs = np.append(cs, c)

            goodf13 = np.append(goodf13, f13)
            goodf123 = np.append(goodf123, f123)
            N = np.append(N, len(cs))
            if len(cs) > 0:
                cmin = np.append(cmin, min(cs))
                cmax = np.append(cmax, max(cs))
            else:
                cmin = np.append(cmin, np.nan)
                cmax = np.append(cmax, np.nan)

    # Save data.
    data = np.vstack((goodf13, goodf123, cmin, cmax, N)).T
    df = pd.DataFrame(data=data, columns=["f13", "f123", "cmin", "cmax", "N"], index=range(1, len(data) + 1))
    df.to_csv(fname)


def plot_compatible_region(fname):
    """ Plots order of magnitude (of free parameter C) over which SsLrpB solutions can be found in the f13-f123 plane,
    the minimal and maximal value of C for compatibility. """

    # Get data
    df = pd.read_csv(fname)

    f13 = df['f13']
    f123 = df['f123']
    cmin = df['cmin']
    cmax = df['cmax']

    # Make f13-f123 grid.
    f13u = np.unique(f13)
    f123u = np.unique(f123)
    f123m, f13m = np.meshgrid(f123u, f13u)

    Nm = np.zeros_like(f13m)
    cminm = np.empty_like(f13m)
    cminm[:] = np.nan
    cmaxm = np.empty_like(f13m)
    cmaxm[:] = np.nan

    # Fill grids of values with data
    for i in range(len(f13m)):
        for j in range(len(f13m[0])):
            cminm[i][j] = cmin[np.logical_and(f13 == f13m[i][j], f123 == f123m[i][j])]
            cmaxm[i][j] = cmax[np.logical_and(f13 == f13m[i][j], f123 == f123m[i][j])]
            Nm[i][j] = np.log10(cmaxm[i][j] / cminm[i][j])

    # Make figure.
    fig = plt.figure(figsize=(10, 3), tight_layout=True)
    ax = fig.add_subplot(131)

    import matplotlib.colors as colors
    ll = ax.pcolor(np.log10(f123m), np.log10(f13m), Nm, vmin=0, vmax=np.max(Nm))
    ll.cmap.set_under('grey')
    fig.colorbar(ll, label=r'Range of solutions [orders of magnitude]')
    ax.set_xlabel(r'$\log(f_{123})$')
    ax.set_ylabel(r'$\log(f_{13})$')

    ct = 1e-3 / 2.4

    Kd1 = 73.5 * ct
    Kd3 = 49.1 * ct
    vmin = 1e-3 * (Kd1 + Kd3)
    vmax = 1e2 * (Kd1 + Kd3)

    ax = fig.add_subplot(132)
    ll = ax.pcolor(f123m, f13m, cminm, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    ll.cmap.set_under('grey')
    fig.colorbar(ll, label=r'$C_{min}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$f_{123}$')
    ax.set_ylabel(r'$f_{13}$')

    ax = fig.add_subplot(133)
    ll = ax.pcolor(f123m, f13m, cmaxm, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    ll.cmap.set_under('grey')
    fig.colorbar(ll, label=r'$C_{max}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$f_{123}$')
    ax.set_ylabel(r'$f_{13}$')

    # Save and show figure.
    # plt.savefig('SsLrpBparams.pdf')
    plt.show()


def plot_bistable_compatible_regions(folder, oscillation_files=0, project=Project.PAPER):
    """ Make figure of bistable region and SsLrpB compatible region."""

    # Generate grid.
    f123 = np.linspace(-3, 2, 51)
    f13 = np.linspace(-3, 2, 51)
    f123, f13 = np.meshgrid(f123, f13)
    T = np.empty_like(f123)
    T[:] = np.nan
    H = np.empty_like(f123)
    H[:] = np.nan
    N = np.zeros_like(f123)

    # Loop over files with solutions and get data.
    for i in range(len(f123)):
        for j in range(len(f123[0])):
            print(f123[i, j], f13[i, j])
            sf123 = int(f123[i][j] * 10)
            if sf123 < 0:
                sf123 = 'm%d' % -sf123
            else:
                sf123 = '%d' % sf123
            sf13 = int(f13[i][j] * 10)
            if sf13 < 0:
                sf13 = 'm%d' % -sf13
            else:
                sf13 = '%d' % sf13
            if os.path.isfile(folder + 'f123-%s-f13-%s.txt' % (sf123, sf13)):
                _, H, _, _, _, Tij, minf, maxf, nonmono = np.loadtxt(folder + 'f123-%s-f13-%s.txt' % (sf123, sf13),
                                                                     unpack=True)
                Tij[H < 50] = np.nan

                T[i][j] = np.nanmax(Tij)
                N[i][j] = sum(np.logical_and(nonmono == 1, np.logical_and(minf < 0.5, maxf > 2.0))) * 100 / len(
                    nonmono)  # in percentage

    # ------------------------------------------------------------------------------------------------------------
    # SEPARATE FIGURES COMPATIBLE AND BISTABLE SOLUTIONS
    # ------------------------------------------------------------------------------------------------------------

    # Create figures with data separated to check for overlap.
    font = {'family': 'Arial', 'size': 10}
    mpl.rc('font', **font)
    fig = plt.figure(figsize=(5.2, 4), tight_layout=True)
    gs = mpl.gridspec.GridSpec(2, 2, height_ratios=[15, 2])  # 1 -> space for legend

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    axcbar1 = plt.subplot(gs[1, 0])
    axcbar2 = plt.subplot(gs[1, 1])

    # Figure with SsLprB compatible solutions.
    N[N == 0] = np.nan
    T = T / 3600  # convert minutes to hours
    ll = ax2.pcolor(f123, f13, N, vmin=0, vmax=100, cmap='Reds')
    ll.cmap.set_under('white')
    fig.colorbar(ll, cax=axcbar2, label=r'Number of solutions', orientation='horizontal')

    # Figure with bistable solutions.
    cmgreen = plt.cm.get_cmap("Greens_r")
    ll = ax1.pcolor(f123, f13, T, cmap=cmgreen, vmin=np.nanmin(T), vmax=np.nanmax(T))
    ll.cmap.set_under('white')
    fig.colorbar(ll, cax=axcbar1, label=r'Induction time [h]', orientation='horizontal')

    # Axes lables, limits and ticks.
    for ax in [ax1, ax2]:
        ax.set_xlabel(r'$\log(f_{123})$', fontsize=10)
        ax.set_ylabel(r'$\log(f_{13})$', fontsize=10, rotation=0)
        ax.yaxis.set_label_coords(-0.02, 1.05)
        ax.set_xlim([-3.5, 2.5])
        ax.set_ylim([-3.5, 2.5])

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)
    #plt.savefig('SsLrpBparamsbistab.pdf')

    # ------------------------------------------------------------------------------------------------------------
    # COMBINED FIGURE
    # ------------------------------------------------------------------------------------------------------------

    # Create the combined figure
    if project == Project.PAPER:
        fig = plt.figure(figsize=(5.2, 3), tight_layout=True)
        gs = mpl.gridspec.GridSpec(6, 2, height_ratios=[1, 4, 1, 4, 4, 1], width_ratios=[3, 2], hspace=0.2)

        ax = plt.subplot(gs[:-1, 0])
        axcbar1 = plt.subplot(gs[0, 1])
        axcbar2 = plt.subplot(gs[2, 1])
        axleg = plt.subplot(gs[4:, 1])
    elif project == Project.SLIDE:
        font = {'family': 'Arial', 'size': 12}

        mpl.rc('font', **font)
        fig = plt.figure(figsize=(8 / 2.54, 7 / 2.54), tight_layout=True)

        ax = plt.subplot(111)
        axcbar1 = 0
        axcbar2 = 0
        axleg = 0

    # Masking to be able to plot both colormaps in one figure, you need to make sure not to hide any found solutions.
    # A check can be made by looking at the separate figures.
    T_ = np.ma.masked_array(T, np.logical_and(f123 <= -0.3, f13 < 1.1))
    N_ = np.ma.masked_array(N, np.logical_and(f123 > -0.3, f13 > 1.1))

    # Plot the SsLrpB compatible solutions.
    ll = ax.pcolor(f123, f13, N, vmin=0, vmax=100, cmap='Reds')
    ll.cmap.set_under('white')
    if axcbar2 != 0:
        fig.colorbar(ll, cax=axcbar2, label='Percentage of Ss-LrpB \ncompatible solutions', orientation='horizontal',
                     ticks=[0, 50, 100, 150, 200])

    # Plot the bistable solutions.
    cmsummer = plt.cm.get_cmap("summer")
    ll = ax.pcolor(f123, f13, T_, cmap=cmsummer, vmin=np.nanmin(T_), vmax=np.nanmax(T_))
    ll.cmap.set_under('white')
    if axcbar1 != 0:
        fig.colorbar(ll, cax=axcbar1, label='Induction time\n of bistable solutions [h]', orientation='horizontal')

    # Axes labels and limits.
    ax.set_xlabel(r'$\log(f_{123})$', horizontalalignment='right', x=1.0)
    ax.set_ylabel(r'$\log(f_{13})$', horizontalalignment='right', x=0, y=0.5)  # , rotation=0)
    # ax.yaxis.set_label_coords(-0.02, 1.05)
    ax.set_xlim([-3.1, 2.1])
    ax.set_ylim([-3.1, 2.1])

    # Add oscillatory solutions.
    if oscillation_files != 0:
        data = pd.read_csv(oscillation_files, index_col=0)

        idx = 347
        f123chosen1 = data.loc[idx]['f123']
        f13chosen1 = data.loc[idx]['f13']
        idx = 1250
        f123chosen2 = data.loc[idx]['f123']
        f13chosen2 = data.loc[idx]['f13']

        cut = np.logical_not(np.isnan(data['period']))

        data = data[cut]
        f123 = np.log10(data['f123'])
        f13 = np.log10(data['f13'])

        # Add oscillating solutions.
        ax.scatter(f123, f13, color='b', marker='o', s=np.ones_like(f123) * 0.1,
                   label='oscillatory solutions, \nSs-LrpB compatible')
        ax.scatter(np.log10(f123chosen1), np.log10(f13chosen1), color='r', marker='*', s=[15],
                   label='selected oscillatory solution A')
        ax.scatter(np.log10(f123chosen2), np.log10(f13chosen2), color='firebrick', marker='*', s=[15],
                   label='selected oscillatory solution B')

        # Legend.
        if axleg != 0:
            handles, labels = ax.get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            leg = axleg.legend(by_label.values(), by_label.keys(), ncol=1, loc=10, frameon=False, handletextpad=0.1)
            axleg.axis('off')
            leg.legendHandles[0]._sizes = [10]
            leg.legendHandles[1]._sizes = [10]

    # Customize tick size.
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(10)  # 18)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(10)  # 18)
    gs.tight_layout(fig)

    # Save and show figure.
    plt.savefig('Figures/SsLrpBparamsbistab-new.eps')
    plt.show()


def main():
    #percentage_SsLrpB_compatible_solutions('bistability/SsLrpB/compatible_solutions.csv')
    #plot_compatible_region('bistability/SsLrpB/compatible_solutions.csv')
    plot_bistable_compatible_regions('bistability/SsLrpB/grid/', 'oscillation/SsLrpB/variables.csv')
    return


if __name__ == "__main__":
    main()
