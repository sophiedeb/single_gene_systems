#!/usr/bin/env python

"""figures_concepts.py: Make the figures of oscillation, bistability, different response curves."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from SGS import *
from matplotlib import gridspec
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum


class Project(Enum):
    PAPER = 1
    SLIDE = 2


def cm2inch(value):
    return value / 2.54


def response_curve(x, phimax, alpha, nA, nR, KA, KR, base, leak):
    """ Return transcription rate with fast dimerization. Curve is generated through generalization of Hill curves."""

    # Increasing Hill curve if activation.
    if nA > 0:
        fact = (alpha * x ** 2) ** nA / (KA ** nA + (alpha * x ** 2) ** nA)
    else:
        fact = 1

    # Decreasing Hill curve if repression.
    if nR > 0:
        frep = 1 / (1 + (alpha * x ** 2 / KR) ** nR)
    else:
        frep = 1

    # Total response is multiplication of the activation and repression.
    return phimax * (base + (1 - base) * fact) * ((1 - leak) * frep + leak)


def concept_figure(project=Project.PAPER):
    """ Save figure with the concept ideas (oscillation and bistability)."""

    # Set font, fontsize, figures and axes depending on the project.
    if project == Project.PAPER:
        font = {'family': 'Arial', 'size': 10}
        mpl.rc('font', **font)
        mpl.rc('legend', handlelength=1)

        fig = plt.figure('Concepts', figsize=(5.5, cm2inch(8)), tight_layout=True)
        gs = gridspec.GridSpec(2, 2, width_ratios=[2, 3], height_ratios=[4, 1], hspace=0.4, wspace=0.3)
        axosc = fig.add_subplot(gs[0])
        axbs = fig.add_subplot(gs[1])
        axoscleg = fig.add_subplot(gs[2])
        axbsleg = fig.add_subplot(gs[3])
    elif project == Project.SLIDE:
        mpl.rc('font', family='Arial', size=16)
        mpl.rc('legend', handlelength=1)

        figbs = plt.figure('Bistability', figsize=(5, 4), tight_layout=True)
        figos = plt.figure('Oscillation', figsize=(5, 4), tight_layout=True)
        axbs = figbs.add_subplot(111)
        axosc = figos.add_subplot(111)
        axbsleg = 0
        axoscleg = 0

    # ------------------------------------------------------------------------------------------------------------
    # FIGURE OSCILLATION
    # ------------------------------------------------------------------------------------------------------------

    # Set A label.
    if project == Project.PAPER:
        wl = -0.15  # width distance label
        hl = 1.1  # height distance label
        fslabels = 10

        axosc.text(wl, hl, 'A', transform=axosc.transAxes, va='top', ha='right', fontsize=fslabels)

    # Set oscillator
    params = {'f12': 0.01036, 'f2': 0.01036, 'vol': 1, 'gammam': 0.8806, 'f1': 100.0, 'DNAtot': 1, 'beta': 2.786,
              'ku2': 0.01, 'kb1': 0.9047, 'taumRNA': 0, 'ku1': 56.28, 'alphadiss': 1.434, 'kb2': 0.002226,
              'gammad': 0.001467, 'bcoop': 1.0, 'ucoop': 0.1737, 'alphaass': 0.9359, 'gammamRNA': 5.232, 'taum': 0,
              'phi0': 2.032}
    ddo = DDS(params)

    # Plot nullclines
    ddo.plot_nullclines(axosc)

    # Calculate local velocities in phase plane.
    dm = []
    dmRNA = []
    for m in np.linspace(0,10,7):
        for mRNA in np.linspace(0,5,6):
            qss = ddo.quasi_steady_state(m)

            var = []

            for x in ddo.ALL_VARS:
                if x == 'mRNA':
                    var += [mRNA]
                else:
                    var += [qss[x]]

            eq = ddo.eqns(var)
            dm += [eq[ddo.ALL_VARS == 'm'][0]]
            dmRNA += [eq[ddo.ALL_VARS == 'mRNA'][0]]

    # Rescale local velocities.
    l = np.sqrt(np.array(dm)**2 + np.array(dmRNA)**2)
    dm /= 1.5*max(l)
    dmRNA /= 1.5*max(l)

    # Plot local velocities as arrows.
    i = 0
    for m in np.linspace(0,10,7):
        for mRNA in np.linspace(0,5,6):
            axosc.arrow(m, mRNA, dm[i], dmRNA[i],
                        head_width=0.25, head_length=0.25, linewidth=0.5, color='grey')
            i += 1

    # Set initial condition for oscillator.
    ss = ddo.steady_state()
    hf = {'m': lambda t: ss['m'][0], 'd': lambda t: 1.3 * ss['d'][0], 'mRNA': lambda t: ss['mRNA'][0],
          'DNA0': lambda t: ss['DNA0'][0], 'DNA1': lambda t: ss['DNA1'][0], 'DNA2': lambda t: ss['DNA2'][0]}

    # Plot time series in phase plane.
    t_f = 1000
    dt = 0.5
    dde = ddo.time_series(t_f, dt, hf=hf)
    ddo.plot_timelapse(0, dde, t_f, dt, DNA=False, axncm=axosc)

    # Set limits, ticks and axis labels.
    axosc.set_xlim([-0.5, 10.5])
    axosc.set_ylim([-0.3, 6])
    axosc.set_xticks([])
    axosc.set_yticks([])
    axosc.set_xlabel('Protein concentration (a.u.)', rotation=0, horizontalalignment='right', x=1.0)
    if project == Project.PAPER:
        axosc.set_ylabel('mRNA concentration (a.u.)')
        axosc.legend_.remove()
        fig.align_labels([axbs, axosc])
    elif project == Project.SLIDE:
        axosc.set_ylabel('mRNA concentration (a.u.)', rotation=0, horizontalalignment='left', y=1.0)

    # Set legend.
    if axoscleg != 0:
        h2, l2 = axosc.get_legend_handles_labels()
        axoscleg.axis('off')
        axoscleg.legend(h2, l2, loc='upper center', frameon=False, labelspacing=0.15)
    else:
        axosc.legend_.remove()

    # ------------------------------------------------------------------------------------------------------------
    # FIGURE BISTABILITY
    # ------------------------------------------------------------------------------------------------------------

    # Set label.
    if project == Project.PAPER:
        axbs.text(wl, hl, 'B', transform=axbs.transAxes, va='top', ha='right', fontsize=fslabels)

    x = np.linspace(0, 800, 10000)

    # Set the parameters.
    alpha = 1  # alphaass/alphadiss
    gammam = 0.02
    KA = 250.0 ** 2
    KR = 400 ** 2
    phimax = 40
    base = 0.01
    leak = 0
    nR = 1
    nA = 2

    # Generate different response functions with no, positive, negative and mixed feedback.
    y = response_curve(x, phimax, alpha, nA, nR, KA, KR, base, leak)  # no feedback
    y2 = response_curve(x, 12, alpha, nA, 0, 180 ** 2, KR, base, leak)  # positive feedback
    y3 = response_curve(x, phimax, alpha, 0, nR, KA, KR, base, leak)  # negative feedback
    y4 = [y2[-1] for _ in range(len(x))]  # mixed feedback

    # Plot the different response curves.
    axbs.plot(x, y4, label='no feedback')
    axbs.plot(x, y2, label='positive feedback')
    axbs.plot(x, y3, label='negative feedback')
    axbs.plot(x, y, label='mixed feedback')
    axbs.plot(x, gammam * x, label='degradation/\ndilution')

    # Vertical dotted lines at L, I, H values.
    axbs.axvline(x=15, color='k', linewidth=0.5, linestyle=':')
    axbs.axvline(x=125, color='k', linewidth=0.5, linestyle=':')
    axbs.axvline(x=605, color='k', linewidth=0.5, linestyle=':')

    # Legend.
    h1, l1 = axbs.get_legend_handles_labels()
    if axbsleg != 0:
        axbsleg.axis('off')
        axbsleg.legend(h1, l1, ncol=2, loc='upper center', frameon=False, labelspacing=0.15)

    # Ticks and axes labels.
    axbs.set_xticks([])
    axbs.set_yticks([])
    axbs.set_xticks([15, 125, 605])
    axbs.set_xticklabels(['$L$', '$I$', '$H$'])
    axbs.set_yticks([0])
    axbs.set_xlabel('Protein concentration (a.u.)', rotation=0, horizontalalignment='right', x=1.0)
    if project == Project.PAPER:
        axbs.set_ylabel('Response function (a.u.)')
    else:
        axbs.set_ylabel('Response function (a.u.)', rotation=0, horizontalalignment='left', y=1.0)

    axbs.set_ylim([0, 40])

    # ------------------------------------------------------------------------------------------------------------
    # SAVE AND SHOW FIGURE(S)
    # ------------------------------------------------------------------------------------------------------------

    if project == Project.PAPER:
        gs.tight_layout(fig)
        plt.savefig('Figures/concepts.eps')
    elif project == Project.SLIDE:
        figos.savefig('Figures/concept-oscillation2.png')
        figbs.savefig('Figures/concept-bistability2.png')

    plt.show()


def nonmonotonic_response_figure():
    """ Generate figure with a nonmonotonic type 1 response curve."""

    # Set font and font size.
    font = {'family': 'Arial', 'size': 10}
    mpl.rc('font', **font)
    mpl.rc('legend', handlelength=1)

    # Set parameters.
    alpha = 1  # alphaass/alphadiss;
    KA = 250.0 ** 2
    KR = 400 ** 2  # 17.0;
    phimax = 40  # 33*4.5 #60.0;
    base = 0.01
    leak = 0
    nR = 1
    nA = 2
    x = np.linspace(0, 800, 10000)

    # Generate response curve.
    y = response_curve(x, phimax, alpha, nA, nR, KA, KR, base, leak)

    # Plot and save response curve.
    plt.figure(figsize=(1.7, 1.5))
    plt.plot(x, y, c='k')
    plt.xticks([])
    plt.yticks([])
    plt.yticks([])
    plt.xlabel('Protein concentration')
    plt.ylabel('Response function')
    plt.title('Non-monotonic response', size=10)
    plt.savefig('Figures/nonmonotonicresponse.pdf')

    plt.show()


def bistability_oscillation_timeseries_figure(project=Project.PAPER):
    """ Generate figure of timeseries of bistability and oscillations."""

    # Set font and figure sizes.
    if project == Project.PAPER:
        # Set font and font size.
        font = {'family': 'Arial', 'size': 10}
        mpl.rc('font', **font)
        mpl.rc('legend', handlelength=1)

        # Set figure size and subplots.
        fig = plt.figure('bistability-oscillation', figsize=(1.7, 3), tight_layout=True)
        gs = gridspec.GridSpec(2, 1, hspace=0.1)

        axbs = fig.add_subplot(gs[0])
        axos = fig.add_subplot(gs[1])
    elif project == Project.SLIDE:
        # Set figure sizes.
        figbs = plt.figure('bistability', figsize=(5, 4), tight_layout=True)
        figos = plt.figure('oscillation', figsize=(5, 4), tight_layout=True)

        axbs = figbs.add_subplot(111)
        axos = figos.add_subplot(111)

    # ----------------------------------------------------------------------
    # BISTABILITY
    # ----------------------------------------------------------------------

    # Define the bistable system.
    params = {'DNAtot': 1, 'taumRNA': 0.0, 'beta': 0.09062, 'f12': 0.06019, 'f23': 0.05841, 'bcoop123': 36.76,
              'ucoop23': 0.08546, 'f1': 0.7795, 'bcoop12': 0.004709, 'bcoop13': 3.622, 'f123': 3.422, 'gammam': 0.1,
              'ku2': 117.5, 'alphaass': 0.01, 'bcoop23': 16.13, 'ucoop12': 0.3368, 'kb3': 1.0, 'f13': 0.003407,
              'alphadiss': 1.0, 'gammamRNA': 0.01, 'f3': 1.956, 'vol': 1, 'phi0': 1.0, 'ucoop123': 0.2904, 'kb2': 1.0,
              'kb1': 1.0, 'ucoop13': 0.6251, 'f2': 0.001173, 'taum': 0.0, 'gammad': 0.001, 'ku1': 1915.0, 'ku3': 0.386}
    ddds = DDDS(params)

    # Do time series starting from different initial conditions.
    starts = [0, 20, 40, 60, 80, 110, 150, 190, 230, 250]
    ss = ddds.steady_state()
    for i in starts:
        # Set the initial conditions.
        if i > 75:
            hf = {'m': lambda t: i, 'd': lambda t: ddds.quasi_steady_state(i)['d'],
                  'mRNA': lambda t: max(0, int(ss['mRNA'][2])), 'DNA0': lambda t: 1, 'DNA1': lambda t: 0,
                  'DNA2': lambda t: 0, 'DNA3': lambda t: 0, 'DNA12': lambda t: 0, 'DNA13': lambda t: 0,
                  'DNA23': lambda t: 0, 'DNA123': lambda t: 0}
        else:
            hf = {'m': lambda t: i, 'd': lambda t: ddds.quasi_steady_state(i)['d'],
                  'mRNA': lambda t: max(0, int(ss['mRNA'][0])), 'DNA0': lambda t: 1, 'DNA1': lambda t: 0,
                  'DNA2': lambda t: 0, 'DNA3': lambda t: 0, 'DNA12': lambda t: 0, 'DNA13': lambda t: 0,
                  'DNA23': lambda t: 0, 'DNA123': lambda t: 0}

        # Generate and plot a time series.
        t_f = 500
        dt = 0.1
        dde = ddds.time_series(t_f, dt, hf=hf)
        ddds.plot_timelapse(axbs, dde, t_f, dt, DNA=False, legend=(i == 0), bw=True, only='m')

    # Ticks, axes labels, limits and legend.
    axbs.set_xticks([])
    axbs.set_yticks([])
    axbs.set_xlabel('Time', x=1, horizontalalignment='right')
    axbs.set_ylabel('Copy number')  # , rotation=0, y=1, horizontalalignment = 'left')
    axbs.set_xlim([0, 150])
    axbs.legend_.remove()

    # Set title/save figure.
    if project == Project.PAPER:
        axbs.set_title('Bistability', size=10, y=0.3, x=0.95, ha='right')
        axbs.axes.get_xaxis().set_visible(False)
    elif project == Project.SLIDE:
        figbs.savefig("Figures/bistability3.pdf")

    # ----------------------------------------------------------------------
    # OSCILLATION
    # ----------------------------------------------------------------------

    # Set the oscillator.
    params = {'f12': 0.01036, 'kb1': 0.9047, 'phi0': 2.032, 'f2': 0.01036, 'gammad': 0.001467, 'taumRNA': 0, 'vol': 1,
              'beta': 2.786, 'DNAtot': 1, 'kb2': 0.002226, 'taum': 0, 'alphadiss': 1.434, 'gammamRNA': 5.232,
              'gammam': 0.8806, 'ucoop': 0.1737, 'alphaass': 0.9359, 'ku1': 56.28, 'ku2': 0.01, 'bcoop': 1.0,
              'f1': 100.0}
    dds = DDS(params)

    # Set the initial conditions.
    hf = {'m': lambda t: 0, 'd': lambda t: 0, 'mRNA': lambda t: 1.0, 'DNA0': lambda t: 1.0, 'DNA1': lambda t: 0,
          'DNA2': lambda t: 0, }

    # Perform a time series and plot it.
    t_f = 5000
    dt = 0.1
    dde = dds.time_series(t_f, dt, hf=hf)
    dds.plot_timelapse(axos, dde, t_f, dt, DNA=False, legend=(i == 0), bw=True, only='m')

    # Ticks, axes labels and limits.
    axos.set_xticks([])
    axos.set_yticks([])
    axos.set_xlabel('Time', x=1, horizontalalignment='right')
    if project == Project.PAPER:
        axos.set_ylabel('Copy number')
    elif project == Project.SLIDE:
        axos.set_ylabel('Copy number', rotation=0, y=1, horizontalalignment='left')
    axos.set_xlim([0, 1000])

    # Set title and save figure.
    if project == Project.PAPER:
        axos.set_title('Oscillation', size=10, y=0.7, x=0.95, ha='right')
        # ax.legend_.remove()
        fig.savefig("Figures/bistabilityoscillation.pdf")
    elif project == Project.SLIDE:
        figos.savefig("Figures/oscillation3.pdf")

    plt.show()


def nonmonotonicities():
    """ Generate figure with different types of nonmonotonicities. """

    # Font, font size and figure size.
    mpl.rc('font', family='Arial', size=10)
    fig = plt.figure(figsize=(5.2, 0.8), tight_layout=True)  # 2,6
    gs = gridspec.GridSpec(1, 5, wspace=0.05)

    # Set parameters.
    alpha = 1  # alphaass/alphadiss;
    KA = 250.0 ** 2
    KR = 400 ** 2
    phimax = 40
    base = 0.01
    leak = 0
    nR = 1
    nA = 2

    # Generate responses of different types.
    x = np.linspace(0, 500, 1000)
    y = response_curve(x, 12, alpha, nA, 0, 180 ** 2, KR, base, leak) # monotonic increasing
    y1 = y[::-1] # monotonic decreasing
    x = np.linspace(0, 1000, 1000)
    y2 = response_curve(x, phimax, alpha, nA, nR, KA, KR, base, leak) # non-monotonic type 1
    x = np.linspace(0, 2500, 1000)
    y3 = -response_curve(x, phimax, alpha, nA, nR, KA, KR, 0.1, leak) # non-monotonic type 2
    x = np.linspace(0, 1500, 1000)
    yh = response_curve(x/1.1, phimax, alpha, nA, nR, KA, KR, 0.2, leak)
    yh2 = [1 - i * 0.004 if i < 100 else 1 - (200 - i) * 0.004 if i < 200 else 1 for i in range(1000)]
    y4 = -yh * yh2 # non-monotonic type 3
    responses = np.vstack((y, y1, y2, y3, y4))

    # Plot all responses.
    titles = ['Increasing', 'Decreasing', 'Type 1', 'Type 2', 'Type 3']
    colors = ['b', 'b', 'r', 'g', 'orange']
    for i in range(len(titles)):
        ax = fig.add_subplot(gs[i])
        # ax.set_title(titles[i)
        ax.plot(x, responses[i], colors[i])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')

    gs.tight_layout(fig)

    # Save and show figure.
    plt.savefig('nonmonotonicities.pdf')
    plt.show()


def main():
    project = Project.PAPER

    concept_figure(project)
    #nonmonotonic_response_figure()
    #bistability_oscillation_timeseries_figure(project)
    #nonmonotonicities()


if __name__ == "__main__":
    main()
