#!/usr/bin/env python

"""SGS.py: Mother class single gene system (SGS) and physiological range."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from scipy.optimize import brentq
from scipy.integrate import odeint
from enum import Enum
import pandas as pd  # for daughter classes


class Polymer(Enum):
    MONOMER = 1
    DIMER = 2


class System(Enum):
    RANDOM = 1
    SSLRPB = 2
    UNREGULATED = 3


def is_int(s):
    "Checks whether string arg is an int."

    try:
        int(s)
        return True
    except ValueError:
        return False


def physiological_range(par):
    # PHI0 CHANGE FOR BISTABILITY, BINDING, UNBINDING

    # Parameter ranges for different binding sites are equivalent.
    equivalences = ['f', 'kb', 'ku', 'K', 'omega']
    for equivalent_par in equivalences:
        if (par.startswith(equivalent_par) and (
                is_int(par[len(equivalent_par):]) or par[len(equivalent_par):] == 'm' or par[
                                                                                         len(equivalent_par):] == 'd')):
            par = equivalent_par
    if (par[1:5] == 'coop' and (
                    (par[0] == 'u' or par[0] == 'b') and len(par) == 5 or is_int(par[5:]) or par[5:] == 'm' or par[
                                                                                                               5:] == 'd')):
        par = 'coop'

    ranges = {'beta': (0.01, 120), 'gammam': (1e-3, 1), 'gammamRNA': (1e-3, 10.0), 'gammad': (5e-4, 0.1),
              'alphaass': (1e-2, 1), 'alphadiss': (1e-3, 1e3), 'phi0': (1e-2, 1e1), 'fa': (1, 100), 'fr': (1e-3, 1),
              'f': (1e-3, 100), 'kb': (0.001, 100), 'ku': (0.01, 1000.0), 'K': (1e-6, 1e4), 'coop': (0.05, 20),
              'Kdim': (0.5, 2), 'omega': (0.05 / 20, 20 / 0.05), 'DNAtot': (0, np.inf), 'taum': (0, np.inf),
              'taumRNA': (0, np.inf)}
    return ranges[par]


class ode23:
    """ Class created as an equivalent to dde23 for delay differential equations."""

    # noinspection PyDefaultArgument
    def __init__(self, oscillator, conc, t, reaccount=[]):
        conc[conc < 0.0] = 0.0  # concentrations cannot be negative, integration errors can become negative
        self.t = t
        self.reaccount = reaccount

        if len(conc[0]) == 2:  # fastdimer
            self.sol = {'mRNA': conc[:, 0], 'm': conc[:, 1]}
        else:
            self.sol = {var: conc[:, i] for i, var in enumerate(oscillator.ALL_VARS)}
            self.sol[oscillator.DNA_STATES[-1]] = oscillator.DNAtot - sum(
                [self.sol[dna] for dna in oscillator.DNA_STATES[:-1]])


class SG:
    def __init__(self):
        return

    # --------------------------------------------------------------------------------------------------------
    # PARAMETER FUNCTIONS
    # --------------------------------------------------------------------------------------------------------

    def set_parameters(self, params):
        """ Set parameters of the system.

        Only parameters of parameterlist are accepted, other parameters are initialized by None.
        Automatic conversion from binding constants to unbinding rates and from cooperativity to unbinding cooperativity.
        The steady state is cleared, because different parameters will lead to a different steady state.
        """

        return

    def get_parameters(self):
        """ Returns dictionary of parameters. """

        params = {}
        for par in self.PARAMETERS:
            params[par] = getattr(self, par)
        return params

    def __copy__(self):
        """ Makes a copy."""

        newone = type(self)()
        newone.set_parameters(self.get_parameters())

        return newone

    def parameter_line(self, sep="\t"):
        """Returns a string of all parameters."""

        string = ""
        for param in self.PARAMETERS:
            string += "%.3E" % getattr(self, param) + sep
        string = string[:-1] + "\n"  # newline instead of tab at end of string
        return string

    def out_of_range(self):
        """Returns a string line with parameters that are out of physiological range or on border."""

        line = ""

        for par in self.parameters:
            if getattr(self, par) < physiological_range(par)[0]:
                line += "%s low; " % par
            elif getattr(self, par) == physiological_range(par)[0]:
                line += "%s low bound; " % par
            elif getattr(self, par) > physiological_range(par)[1]:
                line += "%s high; " % par
            elif getattr(self, par) == physiological_range(par)[1]:
                line += "%s high bound; " % par

        return line

    # set parameters to mutant of parent
    # TODO implement mutant as copy
    def small_mutant(self, parent, constraintf, withdelay=False):
        """ Set parameters to the one of parent with one value slightly changed within the constraints by the constraintfunction.
        """

        not_mutable = ['DNAtot']
        if not withdelay:
            not_mutable += ['taumRNA', 'taum']

        list_params = np.copy(self.PARAMETERS)
        list_params = np.delete(list_params, not_mutable)

        # Choose the parameter that is the be mutated.
        par_mutant = np.random.choice(list_params)

        # Set all parameters to parents parameters except for the mutation which is a value between half and twice the parent value.
        for par in self.PARAMETERS:
            if par == par_mutant:
                mutation = 10 ** (np.random.uniform(np.log10(0.5), np.log10(2)))
                new_value = max(min(getattr(parent, par_mutant) * mutation, constraintf(par_mutant)[1]),
                                constraintf(par_mutant)[0])
                setattr(self, par_mutant, new_value)
            else:
                setattr(self, par, getattr(parent, par))

        # Clear the steady state. New parameters result in a new steady state.
        self.ss.clear()
        return

    def large_mutant(self, parent, constraintf, withdelay=False):
        """ Set parameters to the one of parent with one value randomly changed within the constraints by the constraintfunction.
        """

        not_mutable = ['DNAtot']
        if not withdelay:
            not_mutable += ['taumRNA', 'taum']

        list_params = np.copy(self.PARAMETERS)
        list_params = np.delete(list_params, not_mutable)

        # Choose the parameter that is the be mutated.
        par_mutant = np.random.choice(list_params)

        # Set all parameters to parents parameters except for the mutation which is a random value within the constraints.
        for par in self.PARAMETERS:
            if par == par_mutant:
                new_value = 10 ** np.random.uniform(np.log10(constraintf(par_mutant)[1]),
                                                    np.log10(constraintf(par_mutant)[0]))
                setattr(self, par_mutant, new_value)
            else:
                setattr(self, par, getattr(parent, par))

        # Clear the steady state. New parameters result in a new steady state.
        self.ss.clear()
        return

    # --------------------------------------------------------------------------------------------------------
    # EQUATIONS AND TIME SERIES
    # --------------------------------------------------------------------------------------------------------

    def eqns(self, var, t=0.0, ts=1.0, qs=1.0):
        """
        Ordinary differential equations dictating the single gene system.
        """

        return

    @staticmethod
    def eqns_delay():
        """ Delay differential equations.

        Use of a dictionary of functions for package pydelay.
        """

        return

    def time_series(self, t_f, dt, t_i=0, hf={}):
        """ Return time series of system from t_i to t_f with time step dt starting from initial conditions hf."""

        # Solve ordinary or delay differential equations depending on delay.
        if self.taum == 0 and self.taumRNA == 0:

            # Set initial conditions.
            if len(hf) > 0:
                y0 = [hf[var](0) for var in self.ALL_VARS] #(0)
            else:
                # Start with all DNA empty and no proteins.
                y0 = [self.DNAtot if var == 'DNA0' else 0 for var in self.ALL_VARS]

            # Calculate time series and make ode23 object.
            t = np.arange(t_i, t_f, dt)
            ode = ode23(self, odeint(self.eqns, y0, t), t, [])
            return ode

        else:
            # Time scale and quantity scale
            ts = 1
            qs = 1

            # Delay differential equations need package dde23, works only in Python 2.
            dde = dde23(eqns=self.eqns_delay(), params=self.dimensionless_parameters(ts, qs))
            dde.set_sim_params(tfinal=t_f, dtmax=dt, AbsTol=10 ** -6, RelTol=10 ** -3)
            if len(hf) > 0:
                dde.hist_from_funcs(hf, 50)
            else:
                dde.hist_from_funcs(self.hist_func(qs), 50)
            dde.run()
            return dde

    def stochastic_time_series(self, t_f, dt, t_i=0, hf={}, count_reactions=False):
        """
        Return stochastic time series of system from t_i to t_f with time step dt
        starting from initial conditions hf.

        Only implementation for ODE (no delay)
        """

        if self.taum == 0 and self.taumRNA == 0:

            # Set initial conditions.

            # x = np.zeros((1, len(self.ALL_VARS)))
            if len(hf) != 0:
                x = hf.copy()
            else:
                if not len(self.ss) > 0:
                    self.set_steady_state()
                x = self.integer_state({var: self.ss[var][0] for var in self.ALL_VARS})

            t = t_i
            tt = t_i + dt
            xts = x.copy()

            # Reaction matrix.
            reactions = self.reaction_matrix()

            # Set initial reaction count.
            if count_reactions:
                reaction_count = np.zeros(reactions.shape[1])
            else:
                reaction_count = np.full(reactions.shape[1], np.nan)

            # Do time series.
            while t < t_f:
                # Calculate propensities.
                prop = self.propensities(x)

                # Determine next reaction using the propensities.
                proptot = sum(prop)
                reac = np.random.choice(len(prop), p=prop / proptot)
                Km = np.array([i == reac for i in range(len(prop))])

                # Determine time step tau using the propensities.
                ru = np.random.uniform()
                tau = 1 / proptot * np.log(1 / ru)

                # Add reaction to reaction count.
                if count_reactions:
                    reaction_count[reac] += 1

                # Update time.
                t += tau

                # Save state if time point for saving has passed, define new time point for saving.
                while (
                                t >= tt and tt < t_f - dt / 2):  # save before adapting x (correct for K=1, approximation for K>1 in K-leap method)
                    xts = np.vstack((xts, x))
                    tt += dt

                    # Print time.
                    if tt % 1.0 == 0:
                        print(tt)

                # Update state.
                x += np.asarray(reactions.dot(Km))[0]


            # Save time series as ode23 object.
            ode = ode23(self, xts, np.arange(t_i, t_f, dt), reaction_count)

            return ode

    @staticmethod
    def reaction_matrix():
        """ Returns the reaction matrix.

        Columns are chemical processes. Rows are the change of every variable at given process.
        """

        return

    def propensities(self, var):
        """ Propensities for variables = arg."""

        return

    # --------------------------------------------------------------------------------------------------------
    # STEADY STATE, JACOBIAN & EIGENVALUES
    # --------------------------------------------------------------------------------------------------------

    def set_steady_state(self):
        """
        Calculate and set the steady state of the system.
        """
        return

    def steady_state(self):
        """
        Return steady state of system.
        """

        if not len(self.ss) > 0:
            self.set_steady_state()
        return self.ss

    def quasi_steady_state(self, x, polymer=Polymer.MONOMER):
        """ Returns a dictionary with the quasi steady state for x, the monomer or polymer count.

        steady state --> assuming fast (un)binding and dimerization:
        solution to dm/dt = 0 and dmRNA/dt = 0 under the assumption that dDNAi/dt and dd/dt = 0.
        """

        return

    def jacobian(self, var):
        return

    def eigenvalues(self, ax=0, color='k'):
        if not len(self.ss) > 0:
            self.set_steady_state()
        if max(self.ss['m']) > 1e30:
            return np.array([np.nan])
        for i in range(len(self.ss['m'])):
            ss = [self.ss[var][i] for var in self.ALL_VARS]

            J = self.jacobian(ss)

            # add eigenvalues to matrix (eigenvalues on different rows for different steady states)
            if i == 0:
                eigvals = np.linalg.eigvals(J)
            else:
                eigvals = np.vstack((eigvals, np.linalg.eigvals(J)))

            if ax != 0:
                ax.scatter(eigvals.real, eigvals.imag, c=color, marker='o')

        if ax != 0:
            ax.axvline(0, color='k', linestyle=':')
            ax.set_xlabel("Real")
            ax.set_ylabel("Imaginary")

        return eigvals

    # --------------------------------------------------------------------------------------------------------
    # NULLCLINES
    # --------------------------------------------------------------------------------------------------------

    def nullcline_m(self, x, polymer=Polymer.MONOMER):
        """
        Given x the monomer or dimer concentration,
        the mRNA concentration for which the time derivative of the monomer is zero (dm/dt = 0) is returned.
        """

        qss = self.quasi_steady_state(x, polymer)
        # values of mRNA for which dm/dt are 0 :
        if 'd' in self.ALL_VARS:
            m0 = (2 * self.alphaass * pow(qss['m'], 2) - 2 * self.alphadiss * qss['d'] + self.gammam * qss[
                'm']) / self.beta
        else:
            m0 = self.gammam * qss['m'] / self.beta

        return m0

    def nullcline_mRNA(self, x, polymer=Polymer.MONOMER):
        """
        Given x the monomer or dimer concentration,
        the mRNA concentration for which the time derivative of the mRNA is zero (dmRNA/dt = 0) is returned.
        """

        qss = self.quasi_steady_state(x, polymer)
        # values of mRNA for which dmRNA/dt are 0 :
        mRNA0 = self.phi0 * sum([getattr(self, 'f%s' % dna[3:]) * qss[dna] if dna[3:] != '0' else qss[dna] for dna in
                                 self.DNA_STATES]) / self.gammamRNA
        return mRNA0

    # --------------------------------------------------------------------------------------------------------
    # RATES
    # --------------------------------------------------------------------------------------------------------

    def transcription_rate(self, x, polymer=Polymer.MONOMER):
        """"Returns transcription rate given x, the monomer or dimer."""

        qss = self.quasi_steady_state(x, polymer)

        return self.phi0 * sum(
            [getattr(self, 'f%s' % dna[3:]) * qss[dna] if dna[3:] != '0' else qss[dna] for dna in self.DNA_STATES])

    def translation_rate(self, x, polymer=Polymer.MONOMER):
        """Returns translation rate for x, the monomer or dimer."""

        return self.beta * self.transcription_rate(x, polymer) / self.gammamRNA

    # --------------------------------------------------------------------------------------------------------
    # PLOT FUNCTIONS
    # --------------------------------------------------------------------------------------------------------

    def plot_timelapse(self, ax, de, t_f, dt, t_i=0, DNA=False, axncm=0, axncd=0, legend=True, bw=False, only=0):
        if self.taum == 0 and self.taumRNA == 0:
            t = de.t
            sol = de.sol
        else:  # use pydelay
            if dt < (t_f - t_i) / 5000:  # resample for plot if interval is too small
                dt = (t_f - t_i) / 5000
            sol = de.sample(t_i, t_f, dt)
            t = np.arange(t_i, t_f, dt)

        qs = 1.0  # de.params['qs'];
        ts = 1.0  # de.params['ts'];

        if not DNA:
            m = sol['m']
            mRNA = sol['mRNA']
            if only != 0:
                conc = [sol[only]]
                labels = [only]
                sslabel = [only]
            elif self.alphaass != None:
                if 'd' in sol:
                    d = sol['d']
                else:
                    d = self.quasi_steady_state(m)['d']
                conc = [m, d, mRNA]
                labels = ['Monomer', 'Dimer', 'mRNA']
                sslabel = ['m', 'd', 'mRNA']
            else:
                conc = [m, mRNA]
                labels = ['Monomer', 'mRNA']
                sslabel = ['m', 'mRNA']
        elif 'DNA0' in sol:
            ax.set_ylim([0, self.DNAtot * 1.1])
            conc = [sol[state] for state in self.DNA_STATES[:-1]]
            conc += [self.DNAtot - sum(conc)]
            labels = self.DNA_STATES
            sslabel = labels
        else:
            ax.set_ylim([0, self.DNAtot * 1.1])
            m = sol['m']
            mRNA = sol['mRNA']
            qss = self.quasi_steady_state(m)
            conc = [qss[state] for state in self.DNA_STATES]
            labels = self.DNA_STATES
            sslabel = labels
        if ax != 0:
            for i in range(len(conc)):
                if bw:
                    l = ax.plot(t * ts, conc[i], label=labels[i], color='k')
                else:
                    l = ax.plot(t * ts, conc[i], label=labels[i])
                if not len(self.ss) > 0:
                    self.set_steady_state()
                    print('change in plot')
                if sslabel[i] in self.ss:
                    ls = ['--', ':', '-.']
                    for j in range(len(self.ss[sslabel[i]])):
                        ax.axhline(self.ss[sslabel[i]][j] * qs, color=l[0].get_color(), linestyle=ls[j], linewidth=0.8)
            ax.set_xlabel('Time [min]')  # , fontsize=12);
            if not DNA:
                ax.set_ylabel('Copy number')  # , fontsize=12, rotation=0);
                # ax.yaxis.set_label_coords(-0.02, 1.05)
            # ax.set_ylabel('Concentration [%.2EM]'%qs);
            if legend:
                ax.legend(loc=1, fontsize=11, handlelength=1)
            # ax2 = ax.twinx();
            # mn, mx = ax.get_ylim();
            # ax2.set_ylim(mn*self.vol*qs*navo, mx*self.vol*qs*navo);
            # ax2.set_ylabel('Copy number');
            ax.tick_params(axis='both', which='major', labelsize=11)  # , **calfont)
        if axncm != 0:
            # idx = np.where(np.logical_and(abs(m - m[-1]) < 1e-3, abs(mRNA - mRNA[-1]) < 1e-3))[0][-2] - 1
            # axncm.plot(m[idx:], mRNA[idx:], label='timeseries');

            axncm.plot(m[int(2 * len(m) / 3):], mRNA[int(2 * len(mRNA) / 3):], label='timeseries')
            axncm.set_xlim([0, 1.2 * max(m[int(len(m) / 3):])])
        if axncd != 0:
            # axncd.plot(d[len(d) / 3:] / min(self.ss['d']), mRNA[len(mRNA) / 3:] / min(self.ss['mRNA']));
            axncd.plot(d[int(2 * len(d) / 3):], mRNA[int(2 * len(mRNA) / 3):])
            axncd.set_xlim([0, 1.2 * max(d[int(len(d) / 3):])])
        return ax

    def plot_response_curve(self, ax, polymer=Polymer.MONOMER, maxx=0):
        if not len(self.ss) > 0:
            self.set_steady_state()

        if polymer == Polymer.MONOMER and maxx > 0:
            x = np.arange(0, maxx)
            # m = x
        elif polymer == Polymer.MONOMER:
            x = np.linspace(0, 5 * max(self.ss['m']), 1000)
            # m = x
        elif maxx > 0:
            x = np.arange(0, maxx)
            # m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);
        else:
            x = np.arange(0, 2 * max(self.ss['d']))
            # m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);

        ax.plot(x, self.translation_rate(x, polymer))

    def plot_rate_balance(self, ax, polymer=Polymer.MONOMER, maxx=0):
        if not len(self.ss) > 0:
            self.set_steady_state()

        if polymer == Polymer.MONOMER and maxx > 0:
            x = np.linspace(0, maxx, 1000)
            m = x
        elif polymer == Polymer.MONOMER:
            x = np.linspace(0, 5 * max(self.ss['m']), 1000)
            m = x
        elif maxx > 0:
            x = np.linspace(0, maxx, 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass)
        else:
            x = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass)

        ax.plot(x, self.translation_rate(m, polymer=Polymer.MONOMER), label='translation')
        ax.plot(x, self.gammam * m, label='death')
        if self.alphaass != None:
            # ax.plot(x, self.translationrate(m, monomer=True) + 2 * self.alphadiss * self.alphaass * m ** 2 / (self.alphadiss + self.gammad), label='dimerization')
            # ax.plot(x, self.gammam * m + 2 * self.alphaass * m ** 2, label='dimerization')
            ax.plot(x, 2 * self.gammad * self.alphaass / (self.alphadiss + self.gammad) * m ** 2, label='dimerization')
        ax.plot(x, self.translation_rate(m,
                                         polymer=Polymer.MONOMER) - self.gammam * m - 2 * self.gammad * self.alphaass / (
                    self.alphadiss + self.gammad) * m ** 2, label='sum')

        ax.legend()
        ax.set_ylabel('Rates')
        if polymer == Polymer.MONOMER:
            ax.set_xlabel('Monomer concentration')
            for ssm in self.ss['m']:
                ax.axvline(ssm, color='k', linestyle=':')
        else:
            ax.set_xlabel('Dimer concentration')
            for ssd in self.ss['d']:
                ax.axvline(ssd, color='k', linestyle=':')

    def plot_quasi_steady_state(self, ax, polymer=Polymer.MONOMER, maxx=0):
        """ Plot the quasi steady state of the DNA configurations in ax with monomer or dimer on xaxis."""

        # Set steady state.
        if not len(self.ss) > 0:
            self.set_steady_state()

        # Set monomer and x values.
        if polymer == Polymer.MONOMER and maxx > 0:
            x = np.linspace(0, maxx, 1000)
            m = x
        elif polymer == Polymer.MONOMER:
            x = np.linspace(0, 5 * max(self.ss['m']), 1000)
            m = x
        elif polymer == Polymer.DIMER and maxx > 0:
            x = np.linspace(0, maxx, 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass)
        elif polymer == Polymer.DIMER:
            x = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass)

        # Calculate the quasi steady state.
        qss = self.quasi_steady_state(m)

        # Plot the quasi steady state of all DNA configurations.
        for state in self.DNA_STATES:
            ax.plot(x, qss[state] / self.DNAtot, label=state)

        # Set legend, labels and dotted vertical lines at steady states.
        ax.legend()
        ax.set_ylabel('Ratios at steady-state')
        if polymer == Polymer.MONOMER:
            ax.set_xlabel('Monomer concentration')
            for ssm in self.ss['m']:
                ax.axvline(ssm, color='k', linestyle=':')
        else:
            ax.set_xlabel('Dimer concentration')
            for ssd in self.ss['d']:
                ax.axvline(ssd, color='k', linestyle=':')

    def plot_nullclines(self, ax, polymer=Polymer.MONOMER):
        if not len(self.ss) > 0:
            self.set_steady_state()
        if polymer == Polymer.MONOMER:
            x = np.linspace(0, 20 * max(self.ss['m']), 1000)
        else:
            x = np.linspace(0, 20 * max(self.ss['d']), 1000)

        x0 = self.nullcline_m(x, polymer)
        mRNA0 = self.nullcline_mRNA(x, polymer)

        if polymer == Polymer.MONOMER:
            # x /= min(self.ss['m']);
            x0label = 'dm/dt = 0'
            xlabel = r'Monomer copy number $m$'  # [%.2E M]' % min(self.ss['m']);
        else:
            # x /= min(self.ss['d']);
            x0label = 'dm/dt = 0'
            xlabel = r'Dimer copy number'  # [%.2E M]' % min(self.ss['d']);

        # ax.plot(x, x0 / min(self.ss['mRNA']), label=x0label, color='r');
        # ax.plot(x, mRNA0 / min(self.ss['mRNA']), label='dmRNA/dt = 0', color='orange');
        ax.plot(x, x0, label=x0label, color='r')
        ax.plot(x, mRNA0, label='dmRNA/dt = 0', color='orange')
        ax.legend(loc=2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'mRNA copy number')  # [%.2E M]' % min(self.ss['mRNA']));

    def plot_DNA_time_averaged(self, ax, ode, timestep):
        def helper(X, t, ti, tf):
            if ti < t[0]:
                idxi = 0
            else:
                idxi = np.where((t[:-1] <= ti) & (t[1:] > ti))[0][0]
            if tf >= t[-1]:
                idxf = len(t) - 1
            elif tf < t[0]:
                return 0
            else:
                idxf = np.where((t[:-1] <= tf) & (t[1:] > tf))[0][0] + 1
            th = t[idxi:idxf]
            th[0] = ti
            th[-1] = tf
            Xh = X[idxi:idxf]
            res = np.sum([(te - tb) * nX for tb, te, nX in zip(th[:-1], th[1:], Xh[:-1])])
            return res

        t = ode.t

        newt = np.arange(0, t[-1], timestep)

        DNA = {ode.sol[var] for var in self.DNA_STATES[:-1]}
        DNA[self.DNA_STATES[-1]] = np.full(len(t), self.DNAtot)
        for dna_state in self.DNA_STATES[:-1]:
            DNA[self.DNA_STATES[-1]] -= DNA[dna_state]

        TADNA = {}

        for i in self.DNA_STATES:
            TADNA[i] = np.zeros(len(newt))

        for i in range(1, len(newt)):
            tDNA = [helper(DNA[j], t, newt[i - 1], newt[i]) for j in self.DNA_STATES]
            TADNA[self.DNA_STATES[np.argmax(tDNA)]][i] = 1

        for i in self.DNA_STATES:
            ax.fill_between(newt, 0, TADNA[i], edgecolor="none", label=i, step="post")
        ax.legend()
        ax.set_yticks([])
        ax.tick_params(axis='both', which='major', labelsize=18)  # , **calfont)
        ax.set_xlabel('Time (min)', fontsize=12)
        ax.legend(loc='upper center', markerscale=0.35, columnspacing=1, bbox_to_anchor=(0.5, -1.5), ncol=4,
                  fontsize=12)

    # --------------------------------------------------------------------------------------------------------
    # MONOTONICITY
    # --------------------------------------------------------------------------------------------------------

    """
    def monotonic(self):
        if (not len(self.ss) > 0):
            self.set_steady_state();
        x = np.linspace(min(self.steady_state()['m']), 5 * max(self.steady_state()['m']), 2000);
        tx = self.translationrate(x);
        if ((max(tx) - tx[-1]) / max(tx) > 0.05):
            return False
        else:
            return True
    """

    def monotonicity_type(self):
        """
        Returns the monotonicity type:

        0 -> monotonic
        1 -> non-monotonic with one maximum
        2 -> non-monotonic with one minimum
        3 -> non-monotonic with a maximum (maxima) and a minimum (minima)
        """

        # Make response curve.
        xx = np.linspace(0, 1000, 1000)
        xx[0] = 0.01
        y = self.transcription_rate(xx, polymer=Polymer.DIMER)

        # Look where response curve has local maximum/minimum (derivative of response curve changes sign).
        dy = y[1:] - y[:-1]
        zero_crossings = np.where(np.diff(np.sign(dy)))[0]

        # Count local extrema and type of extrema to determine the monotonicity type.
        if len(zero_crossings) == 1:
            if dy[zero_crossings[0]] > 0:
                nonmono = 1
            elif dy[zero_crossings[0]] < 0:
                nonmono = 2
        elif len(zero_crossings) > 1:
            nonmono = 3
        else:
            nonmono = 0

        return nonmono

    """
    def nonmonotonicity(self):
        if (not len(self.ss) > 0):
            self.set_steady_state();
        x = np.linspace(min(self.ss['m']), max(self.ss['m']), 100)
        maxtrans = max(self.transcription_rate(x))
        endtrans = self.transcription_rate(x[-1])
        # endtrans = self.transcriptionrate(max(self.ss['m']))
        # if(maxtrans-self.transcriptionrate(x[-1])>0):
        #    return 1 + (maxtrans-endtrans)/maxtrans
        # else:
        return (maxtrans - endtrans) / maxtrans
    """

    # --------------------------------------------------------------------------------------------------------
    # BISTABILITY
    # --------------------------------------------------------------------------------------------------------

    def is_bistable(self):
        """ Returns if system is bistable."""
        # Set steady state.
        if not len(self.ss) > 0:
            self.set_steady_state()
        if len(self.ss['m']) == 3:
            eig = self.eigenvalues()
            stability = np.all(eig.real < 0, axis=1)

            if sum(stability) == 2: #2 stable steady states
                return True
            else:
                print("stab", stability, self.get_parameters())
        return False

    def is_multistable(self):
        """ Returns if system is multistable. Can have oscillations around one or multiple of fixed points."""

        # Set steady state.
        if not len(self.ss) > 0:
            self.set_steady_state()
        return len(self.ss['m']) > 1

    def stochastic_induction_time(self):  # TODO check initial conditions
        """
        Return time it takes for a stochastic time series to go from the intermediate steady state to the high steady state.

        Only implementation for ODE (no delay)
        """

        if self.taum == 0 and self.taumRNA == 0:
            # Set steady state.
            if not len(self.ss) > 0:
                self.set_steady_state()

            # initial condition -> integer version of a perturbation on the deterministic intermediate steady state
            x_det = {var: (self.ss[var] if var in self.DNA_STATES else 1.2 * self.ss[var]) for var in self.ALL_VARS}
            x = self.integer_state(x_det)

            # x = np.array([[1, 0, 0, 0, 0, 0, 0, int(self.ss['mRNA'][1]) + 1, int(self.ss['m'][1]) + 1,
            #               int(self.ss['d'][1]) + 15]])

            t = 0
            tt = 50
            maxt = 5000

            # Reaction matrix.
            reactions = self.reaction_matrix()

            # End state.
            maxd = 0.95 * max(self.steady_state()['d'])

            xs = np.copy(x)

            # Perform time series until end state is reached or time is larger than max time.
            while t < maxt and x[-1] < maxd:
                # Calculate propensities
                prop = self.propensities(x[0])

                # Determine next reaction and time step using the propensities.
                proptot = sum(prop)
                reac = np.random.choice(len(prop), p=prop / proptot)
                Km = np.array([i == reac for i in range(len(prop))])
                ru = np.random.uniform()
                tau = 1 / proptot * np.log(1 / ru)

                # Adapt time.
                t += tau

                # Print intermediate time steps.
                if t > tt:
                    print(
                        "time is %d, dimer count is %d (I = %d, H = %d)" % (t, x[-1], self.ss['d'][1], self.ss['d'][2]))
                    tt = 50 * (int(t / 50) + 1)

                # Adapt state.
                x += reactions.dot(Km)

            return t

    """
    def score_bistab(self):
        if (not len(self.ss) > 0):
            self.set_steady_state();

        # score for bistability on 10:
        if (len(self.ss['m']) < 3):
            x = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);

            x = self.translationrate(m, monomer=True) - self.gammam * m - 2 * self.gammad * self.alphaass / (
            self.alphadiss + self.gammad) * m ** 2

            lm = np.r_[False, x[1:] < x[:-1]] & np.r_[x[:-1] < x[1:], False]  # local minima
            if (sum(lm) == 0):
                score = 0;
            elif (sum(lm) > 1):
                print("sum bigger than 1")
                score = 0
            else:
                score = 10.0 / (1 + np.exp(-2 * (x[lm][0] - 10)))
        else:
            scorethreshold = np.tanh((self.ss['d'][1] - self.ss['d'][0]) / 20.0) * np.tanh(
                (self.ss['d'][2] - self.ss['d'][1]) / 20.0)
            scoremax = np.tanh(self.ss['d'][2] / 100.0) * 1 / (1 + np.exp(0.01 * (self.ss['d'][2] - 1500.0)))

            scoremRNAmax = 1 / (1 + np.exp((self.ss['mRNA'][2] - 20.0)))

            d = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * d / self.alphaass);

            x = self.translationrate(m, monomer=True) - self.gammam * m - 2 * self.gammad * self.alphaass / (
                self.alphadiss + self.gammad) * m ** 2;

            posx = x[np.logical_and(d > self.ss['d'][1], d < self.ss['d'][2])]

            scoretl = np.tanh(sum(posx) / len(posx) / 5.0)

            # bigx = x[np.logical_and(d > self.ss['d'][2], d < 1.15*self.ss['d'][2])]
            # scoretl2 = np.tanh(-sum(bigx)/len(bigx)/5.0)

            # _,_,_,_,it = self.parambistab()

            # scoreit = 1/(1+np.exp(0.01*(it-500.0)))

            # lm = np.r_[False, x[1:] > x[:-1]] & np.r_[x[:-1] > x[1:], False]  & np.r_[d > self.ss['d'][1]] & np.r_[d < self.ss['d'][2]] # local maximum between unstable and stable steady state
            # if(len(x[lm])==1):
            #    scoretl = np.tanh(x[lm][0]/5.0);
            # else:
            #    scoretl = 1
            #    print("multiple maxima")
            score = 10 + 2.5 * scorethreshold + 2.5 * scoremax + 2.5 * scoremRNAmax + 2.5 * scoretl  # + 2.5*scoretl + 2.5*scoretl2;
        return 20 - score;
    """

    # TODO: check
    def deterministic_induction_time(self,
                                     stoch=False):
        """ """
        if not len(self.ss) > 0:
            self.set_steady_state()

        if self.is_multistable() and len(self.ss['d']) > 2:

            t_f = 10000
            dt = 1E-2

            hf = {}

            def set_lambda(
                    par):  # need this function otherwise problems with loops over dictionaries of lambda functions
                if par == "d":
                    return lambda t: 1.1 * self.steady_state()[par][1]
                else:
                    return lambda t: self.steady_state()[par][1]

            for par in self.ALL_VARS:
                hf[par] = set_lambda(par)

            ts = self.time_series(t_f, dt, hf=hf)
            tsd = ts.sol['d']

            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # self.plottimelapse(ax, ts, t_f, dt)
            # plt.show()

            it = len(tsd[tsd < 0.90 * self.steady_state()['d'][2]]) * dt
            if False:
                for par in self.ALL_VARS:
                    print(par, hf[par](0))
                print(tsd)
                fig = plt.figure("ts %.3E" % it)
                ax = fig.add_subplot(1, 1, 1)
                self.plot_timelapse(ax, ts, t_f, dt)
                ax.axhline(0.95 * self.steady_state()['d'][2])
                plt.show()

            if stoch:
                itstoch = np.zeros(5)
                for i in range(5):
                    itstoch[i] = self.stochastic_induction_time()

                return it, itstoch[0], itstoch[1], itstoch[2], itstoch[3], itstoch[4]
            else:
                return it
        else:
            return np.nan

    # --------------------------------------------------------------------------------------------------------
    # OSCILLATION
    # --------------------------------------------------------------------------------------------------------

    def oscillation_parameters(self, vars):  # eigenvalues have to have positive real part, otherwise nan
        """ Returns the period + minimum and maximum of the amplitude of all variables in the arg. """

        # Count number of variables
        if isinstance(vars, str):  # Python 3: str, python 2 : basestring
            N = 1
        else:
            N = len(vars)

        # Set steady state if not yet done.
        if not len(self.ss) > 0:
            self.set_steady_state()

        # Return NaN if system is multistable.
        if self.is_multistable():
            return [np.nan] * (1 + 2 * N)

        # Return NaN if system is not oscillatory.
        #if not self.is_oscillatory():
        #    return [np.nan] * (1 + 2 * N)

        # Choose variable with highest steady state for finding period.
        hvar = max(['mRNA', 'm', 'd'], key=lambda x: self.ss[x][0])

        # Starting values for the time series.
        t_f = 3000
        dt = 0.1
        hf = {}

        # Run time series ever longer and with smaller time steps until period can be abstracted.
        while t_f < 1e10:
            # Do time series.
            de = self.time_series(t_f, dt, hf=hf)
            if self.taum == 0 and self.taumRNA == 0:
                sol = de.sol
            else:
                sol = de.sample(0, t_f, dt)

            # Quantity scale qs and time scale ts
            qs = 1
            ts = 1

            # Get locations of the maxima and minima.
            locmax = argrelextrema(sol[hvar], np.greater)[0]
            locmin = argrelextrema(sol[hvar], np.less)[0]

            # To calculate period, we want at least 5 maxima and 5 minima, and the height of these extrema must be
            # the same (less than 1% difference) for regular oscillations.
            # Do not use last elements of locmax/locmin ([-1]), they could be the last elements of the time series
            # and not a real local maximum/minimum of the oscillation.
            if (len(locmax) > 5 and len(locmin) > 5 and (
                        abs(sol[hvar][locmax[-4]] - sol[hvar][locmax[-3]]) < 0.01 * sol[hvar][locmax[-3]]) and (
                        abs(sol[hvar][locmin[-4]] - sol[hvar][locmin[-3]]) < 0.01 * sol[hvar][locmin[-3]]) and (
                        (2 * locmax[-3] - locmax[-2] - locmax[-4]) < 3)):

                def period_ampmin_ampmax(var):
                    # Get the local maxima and minima locations of the variable.
                    if var != hvar:
                        locmax_var = argrelextrema(sol[var], np.greater)[0]
                        locmin_var = argrelextrema(sol[var], np.less)[0]
                    else:
                        locmax_var = locmax
                        locmin_var = locmin

                    # If no good oscillatory behavior for variable, return NaN.
                    if len(locmax_var) < 3 or len(locmin_var) < 3:
                        return [np.nan] * 3

                    # Calculate period, maximum and minimum values from the locations of the extrema.
                    period = dt * ts * (locmax_var[-2] - locmax_var[-3])
                    ampmin = qs * sol[var][locmin_var[-2]]
                    ampmax = qs * sol[var][locmax_var[-2]]

                    return period, ampmin, ampmax

                # Calculate the period, maximum and minimum amplitude of the variable(s).
                if N == 1:  # Python 3: str, python 2 : basestring
                    period, ampmin, ampmax = period_ampmin_ampmax(vars)

                    return [period, ampmin, ampmax]
                else:
                    ospars = []
                    for i, var in enumerate(vars):
                        period, ampmin, ampmax = period_ampmin_ampmax(var)

                        # Add minimal and maximal amplitude to oscillation parameters.
                        if i == 0:
                            period = dt * ts * (locmax[-2] - locmax[-3])
                            ospars += [period, ampmin, ampmax]
                        else:
                            ospars += [ampmin, ampmax]

                    return ospars
            else:
                # If period not found, double length of time series and half the time step, do again until period found
                # or time step becomes too computationally expensive (t_f / dt < 1e6).
                t_f = 2 * t_f
                dt = dt / 2
                hf = {var : (lambda t : sol[var][-1]) for var in self.ALL_VARS}

                if (t_f/dt > 1e6):
                    dt = t_f/1e6

        return [np.nan] * (1 + 2 * N)

    def is_oscillatory(self):
        """ Returns whether the system has Hopf oscillations (can also be oscillations + steady state).

        Checks if any of the eigenvalues has a positive real part and non-zero imaginary part.
        """

        # Set steady state if not yet done.
        if not len(self.ss) > 0:
            self.set_steady_state()

        if len(self.ss['m']) > 1:
            return False

        eig = (self.eigenvalues()).flatten()
        return np.any(eig.real > 0) # np.logical_and(eig.imag != 0,eig.real > 0))

    # --------------------------------------------------------------------------------------------------------
    # AUXILIARY FUNCTIONS
    # --------------------------------------------------------------------------------------------------------

    def integer_state(self, state):
        """ Returns state with only integer numbers for stochastic simulations."""

        int_state = {}
        # Set all variables to the rounded value of the real steady state.
        for var in self.ALL_VARS:
            int_state[var] = np.round(state[var])
        int_state[self.DNA_STATES[-1]] = np.round(self.DNAtot - sum([int_state[dna] for dna in self.DNA_STATES[:-1]]))

        # Check whether the total DNA is right. Otherwise fix it.
        diff_dna = self.DNAtot - sum([int_state[dna] for dna in self.DNA_STATES])
        if diff_dna > 0:
            max_dna = self.DNA_STATES[np.argmax([state[dna] for dna in self.DNA_STATES])]
            int_state[max_dna] += diff_dna

        return int_state

    def hist_func(self, qs=1.0, fastdimer=False):
        """ Default history function, a slight perturbation of the steady state.

        Scaled with quantity scale qs.
        """

        # Set steady state if not yet done.
        if not len(self.ss) > 0:
            self.set_steady_state()

        # Set all DNA states equal to steady states and proteins at 120% of steady state.
        hf = {}
        # if(fastdimer):
        #    hf = {'m': lambda t: self.ss['m'] / qs * 1.1, 'mRNA': lambda t: self.ss['mRNA'] / qs * 0.9}
        # else:
        # Set all DNA states equal to steady states and proteins at 120% of steady state.
        for var in self.ALL_VARS:
            if var not in self.DNA_STATES:
                hf[var] = lambda t: self.ss[var] / qs * 1.2
            else:
                hf[var] = lambda t: self.ss[var] / qs * 1.2
        return hf

    def dimensionless_parameters(self, ts, qs):
        """Returns dimensionless parameter list using timescale ts and quantity scale qs."""

        return

    def extrema(self, t_f, dt, vars, hf={}):
        """
        Returns extrema of timeseries of lenght t_f with time step dt, using variables vars
        starting at initial conditions hf.
        """

        if len(hf) == 0:
            hf = self.hist_func()
        ode = self.time_series(t_f, dt, hf=hf).sol
        # function to remove adjacent duplicates
        remdup = lambda x: np.array([n for i, n in enumerate(x) if i == 0 or n != x[i - 1]])
        ext = {}
        for var in vars:
            if var in self.ALL_VARS:
                # idxvar = self.allvars.index(var);
                odevar = ode[var]
            elif var in self.DNA_STATES:
                odevar = self.DNAtot - sum(
                    [ode[i] for i in self.DNA_STATES[:-1]])  # ode[:-1, 0] - ode[:-1, 1] - ode[:-1, 2];
            else:
                print('variable %s not found' % var)
                return
            odevar = remdup(odevar)
            ode0 = odevar[:-1]
            ode1 = odevar[1:]
            pos = [(np.r_[False, ode1 < ode0] & np.r_[ode0 < ode1, False]) | (
                np.r_[False, ode1 > ode0] & np.r_[ode0 > ode1, False])]  # first and last elements always false
            ext[var] = ode0[pos[0][:-1]]
        return ext

    def scores(self, t_f, dt, hf={}):
        scores = [20 for _ in range(len(self.ALL_VARS) + 1)]
        if False:
            for i in range(len(self.DNA_STATES)):
                scores[i] = self.scorefunction(t_f, dt, self.DNA_STATES[i], hf)
            for i in range(len(self.ALL_VARS) - (len(self.DNA_STATES) - 1)):
                scores[len(self.DNA_STATES) + i] = self.scorefunction(t_f, dt,
                                                                      self.ALL_VARS[len(self.DNA_STATES) - 1 + i], hf)
        else:
            vars = self.DNA_STATES
            for i in range(len(self.DNA_STATES) - 1, len(self.ALL_VARS)):
                vars = np.append(vars, self.ALL_VARS[i])
            scores = self.scorefunction(t_f, dt, vars, hf)
        return scores

    def scoreline(self, t_f, dt, hf={}):
        scores = self.scores(t_f, dt, hf)
        line = ""
        for score in scores:
            line += "%.3E\t" % score
        line += "\n"
        return line

    def scorefunction(self, t_f, dt, var, hf={}):
        if not len(self.ss) > 0:
            self.set_steady_state()
        if type(var) == str:
            ext = self.extrema(t_f, dt, [var], hf)[var]
            S = 20
            S -= 1 / (1 + np.exp((self.ss['mRNA'][0] - 20.0)))
            S -= np.tanh(self.ss['d'][0] / 20.0)
            S -= np.tanh((self.ss['d'][0] - self.ss['m'][0]) / 20.0)
            S -= np.tanh((self.ss['d'][0] - self.ss['mRNA'][0]) / 20.0)
            if len(ext) >= 2:
                S -= 0.5
            if len(ext) >= 5:
                S -= 0.5
            if len(ext) >= 7:
                S -= 0.5
            if len(ext) >= 9:
                S -= 0.5
            if len(ext) > 2:
                if len(ext) <= 11:
                    begin = 1
                else:
                    begin = 1  # len(ext) - 11;
                for i in range(begin, min(11, len(ext) - 1)):  # len(ext) - 1):
                    S = S - abs(ext[i] - ext[i + 1]) / (ext[i] + ext[i + 1]) * min(1., abs(ext[i] - ext[i + 1]))
        else:
            S = [20 for _ in range(len(var))]
            extrema = self.extrema(t_f, dt, var, hf)
            for k in range(len(var)):
                ext = extrema[var[k]]
                S[k] -= 1 / (1 + np.exp((self.ss['mRNA'][0] - 20.0)))
                # if(isinstance(self, MDS)):
                #    S[k] -= np.tanh(self.ss['d'][0] / 20.0);
                #    S[k] -= np.tanh(self.ss['m'][0] / 20.0);
                # else:
                S[k] -= np.tanh(self.ss['d'][0] / 50.0)
                S[k] -= np.tanh((self.ss['d'][0] - self.ss['m'][0]) / 20.0)
                if len(ext) >= 3:
                    S[k] -= 0.5
                if len(ext) >= 5:
                    S[k] -= 0.5
                if len(ext) >= 7:
                    S[k] -= 0.5
                if len(ext) >= 9:
                    S[k] -= 0.5
                if len(ext) > 2:
                    if len(ext) <= 11:
                        begin = 1
                    else:
                        begin = 1  # len(ext) - 10;
                    for i in range(begin, min(11, len(ext) - 1)):
                        S[k] = S[k] - abs(ext[i] - ext[i + 1]) / (ext[i] + ext[i + 1]) * min(1.,
                                                                                             abs(ext[i] - ext[i + 1]))
                        if ext[i] + ext[i + 1] == 0:
                            S[k] = 20.0
                            # S = np.mean(S)
        return S

    @staticmethod
    def _find_x_steady_state(func):
        """ Find steady states of x between 1e-25 and 1e25. """

        # Check whether steady state is not too high or low, else determine interval around steady state(s)
        # and determine steady states by looking for zeroes with function brentq.
        if func(1e25) > 0:
            y = np.array([np.inf])
        elif func(1e-25) < 0:
            y = np.array([0])
        else:
            N = 5000
            x = np.logspace(-25, 25, N)
            zeroes = np.where(np.diff(np.sign(func(x))))[0]

            # If multiple zeroes, there is multistability, all steady states need to be calculated.
            # Otherwise single steady state is calculated.
            if len(zeroes) > 1:
                zeroes = np.append(np.append([0], zeroes), [N - 1])
                y = np.array([brentq(func, (x[zeroes[i]] + x[zeroes[i - 1] + 1]) / 2,
                                     (x[zeroes[i + 1]] + x[zeroes[i] + 1]) / 2, xtol=1e-10, rtol=1e-5, maxiter=200) for
                              i in range(1, len(zeroes) - 1)])
            else:
                y = np.array([brentq(func, 1e-25, 1e25, xtol=1e-10, rtol=1e-6, maxiter=200)])

        return y
