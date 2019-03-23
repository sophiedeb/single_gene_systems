#!/usr/bin/env python

"""SGS.py: Import of all SGS (MDS, 2DS, 3DS), auxiliary functions to generate parameters, lists, path"""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from MDS import *
from DDS import *
from DDDS import *
import time

def random_params(sgs, constraintf, system=System.RANDOM, seed=0):
    """ Return random set of parameters withing the constraints. """

    params = {'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'vol': 1}

    if sgs == MDS:
        if system == System.UNREGULATED:
            params['kbm'] = 0
            params['kbd'] = 0
            params['kum'] = np.inf
            params['kud'] = np.inf
            params['fm'] = 0
            params['fd'] = 0
            par_list = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0']
        elif system != System.RANDOM:
            print("Not a known system for this single gene system.")
            return
        par_list = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'kbm', 'kbd', 'kum', 'kud',
                'fm', 'fd']
    elif sgs == DDS:
        if system != System.RANDOM:
            print("Not a known system for this single gene system.")
            return
        par_list = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'kb1', 'kb2', 'ku1', 'ku2', 'ucoop', 'bcoop', 'f1', 'f2', 'f12']
    elif sgs == DDDS:
        if system == System.SSLRPB:
            ct = 1e-3 / 2.4

            params['Kd1'] = 73.5 * ct
            params['Kd2'] = 0.7 * ct
            params['Kd3'] = 49.1 * ct
            omegas = {12: 4.1, 23: 2.1, 13: 7.9, 123: 3.1}

            par_list = np.array(
                ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'kb1', 'kb2', 'kb3',
                 'bcoop12', 'bcoop23', 'bcoop13', 'bcoop123', 'f1', 'f2', 'f3', 'f12', 'f13', 'f23', 'f123'])

            # shuffle parameters in random order
            np.random.seed(int(time.time()))
            np.random.shuffle(par_list)

            if seed != 0:
                np.random.seed(seed)

            for par in par_list:
                if par.startswith('kb'):
                    idx = int(par[-1])
                    params['kb%d' % idx] = 10 ** np.random.uniform(np.log10(max(constraintf('kb')[0], constraintf('ku')[0] * params['Kd%d' % idx])),
                                                                   np.log10(min(constraintf('kb')[1], constraintf('ku')[1] * params['Kd%d' % idx])))
                elif par.startswith('bcoop'):
                    idx = int(par[5:])
                    params['bcoop%d' % idx] = 10 ** np.random.uniform(np.log10(constraintf('coop')[0] * omegas[idx]), np.log10(constraintf('coop')[1]))
                    params['ucoop%d' % idx] = params['bcoop%d' % idx] / omegas[idx]
                else:
                    params[par] = 10 ** (
                    np.random.uniform(np.log10(constraintf(par)[0]), np.log10(constraintf(par)[1])))

            return params
        elif system != System.RANDOM:
            print("Not a known system for this single gene system.")
            return
        else:
            par_list = np.array(
                ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'kb1', 'kb2', 'kb3', 'ku1',
                 'ku2', 'ku3', 'bcoop12', 'bcoop23', 'bcoop13', 'bcoop123', 'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123',
                 'f1', 'f2', 'f3', 'f12', 'f13', 'f23', 'f123'])

    # Shuffle parameters in random order.
    np.random.seed(int(time.time()))
    np.random.shuffle(par_list)

    # Choose random parameters.
    if seed != 0:
        np.random.seed(seed)

    for par in par_list:
        params[par] = 10 ** (
        np.random.uniform(np.log10(constraintf(par)[0]), np.log10(constraintf(par)[1])))

    return params

def parameter_list_latex(sgs):
    if sgs == MDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'kbm', 'kbd', 'Km',
               'Kd']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_m$', r'$f_d$', r'$k_{bm}$',
                    r'$k_{bd}$', r'$K_{m}$', r'$K_{d}$']
    elif sgs == DDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12', 'kb1',
               'kb2', 'ku1', 'ku2', 'bcoop', 'ucoop']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_1$', r'$f_2$', r'$f_{12}$',
                    r'$k_{b1}$', r'$k_{b2}$', r'$k_{u1}$', r'$k_{u2}$', r'coop$_{b12}$', r'coop$_{u12}$']
    elif sgs == DDDS:
        var = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12', 'f13',
               'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop23', 'bcoop13', 'bcoop123',
               'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']
        varlabel = [r'$\beta$', r'$\gamma_\mathrm{m}$', r'$\gamma_\mathrm{mRNA}$', r'$\gamma_\mathrm{d}$',
                    r'$\alpha_\mathrm{ass}$', r'$\alpha_\mathrm{diss}$', r'$\phi_0$', r'$f_1$', r'$f_2$', r'$f_3$',
                    r'$f_{12}$', r'$f_{13}$', r'$f_{23}$', r'$f_{123}$', r'$k_{b1}$', r'$k_{b2}$', r'$k_{b3}$',
                    r'$k_{u1}$', r'$k_{u2}$', r'$k_{u3}$', r'coop$_{b12}$', r'coop$_{b23}$', r'coop$_{b13}$',
                    r'coop$_{b123}$', r'coop$_{u12}$', r'coop$_{u23}$', r'coop$_{u13}$', r'coop$_{u123}$']
    return var, varlabel

class Concept(Enum):
    OSCILLATION = 1
    BISTABILITY = 2

def foldername(sgs, system, concept):
    """ Return pathname given the concept and system."""

    print(concept == Concept.OSCILLATION)

    if concept == Concept.OSCILLATION:
        print("in os")
        folder = 'oscillation/'
    elif concept == Concept.BISTABILITY:
        folder = 'bistability/'

    if sgs == DDDS and system == System.RANDOM:
        folder += '3DS/'
    elif sgs == DDS and system == System.RANDOM:
        folder += '2DS/'
    elif sgs == MDS and system == System.RANDOM:
        folder += 'MDS/'
    elif sgs == DDDS and system == System.SSLRPB:
        folder += 'SsLrpB/'

    return folder