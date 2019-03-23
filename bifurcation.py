"""bifurcation.py: Bifurcation analysis for oscillating/bistable solutions."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from SGS import *
import multiprocessing
import pandas as pd
import os as pyos
from functools import partial

num_cores = 2

def one_core(idx, sgs = DDDS, folder = '', concept = Concept.OSCILLATION):
    """ For one oscillating/bistable solution, save the bifurcation points for every parameter
     and five points in between."""

    bifpoints_name = folder + "bifurcations/points-%d.csv" % idx
    bif_name = folder + "bifurcations/bifurcation-%d.csv" % idx

    if not pyos.path.exists(bif_name):
        filename = folder + "variables.csv"
        data = pd.read_csv(filename, index_col=0)

        # -------------------------------------------------------------------------------
        # parameters
        # -------------------------------------------------------------------------------

        params = data.loc[idx].to_dict()
        os = sgs(params)

        bifurcation_endpoints(os, bifpoints_name, concept)

        bif_endpoints = pd.read_csv(bifpoints_name, index_col=0)

        bifurcation_sampling(os, bif_name, bif_endpoints, concept)

    return

def bifurcation_endpoints(sgs, fname, concept = Concept.OSCILLATION):
    """ Save left and right bounds of oscillatory/bistable range for all system parameters in csv file.

    Bounds are limited by physiological range.
    """

    # -------------------------------------------------------------------------------
    # bifurcation parameters
    # -------------------------------------------------------------------------------

    # List all bifurcation parameters.
    if isinstance(sgs, DDDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3',
                   'f12', 'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3',
                   'bcoop12', 'bcoop23', 'bcoop13', 'bcoop123', 'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']
    elif isinstance(sgs, DDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12',
               'kb1', 'kb2', 'ku1', 'ku2', 'bcoop', 'ucoop']
    elif isinstance(sgs, MDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm',
               'kbd']

    # -------------------------------------------------------------------------------
    # bifurcation for every parameter
    # -------------------------------------------------------------------------------

    params = sgs.get_parameters()

    bounds = np.zeros((len(parnames),3))

    sgs_run = type(sgs)()

    for i in range(len(parnames)):

        # Set parameters of "running" oscillator to copy of original parameters.
        params_run = { par : params[par] for par in parnames }
        for par in ['taum', 'taumRNA', 'DNAtot']:
            params_run[par] = params[par]

        # -------------------------------------------------------------------------------
        # find bifurcation points, Newton method
        # -------------------------------------------------------------------------------

        # Left bound will be between left side of physiological range and value of orignal oscillator.
        left_bound1 = np.log10(physiological_range(parnames[i])[0])
        left_bound2 = np.log10(params[parnames[i]])

        # If left bound 1 is oscillating, it is the left bound of the range.
        params_run[parnames[i]] = 10**left_bound1
        sgs_run.set_parameters(params_run)
        if ((concept == Concept.OSCILLATION and sgs_run.is_oscillatory()) or
            (concept == Concept.BISTABILITY and sgs_run.is_bistable())):
            left_bound = left_bound1
        else:
            while left_bound2-left_bound1 > 0.01:
                mid = (left_bound2+left_bound1)/2
                params_run[parnames[i]] = 10**mid
                sgs_run.set_parameters(params_run)
                if ((concept == Concept.OSCILLATION and sgs_run.is_oscillatory()) or (
                        concept == Concept.BISTABILITY and sgs_run.is_bistable())):
                    left_bound2 = mid
                else:
                    left_bound1 = mid
            left_bound = (left_bound2+left_bound1)/2

        right_bound1 = np.log10(params[parnames[i]])
        right_bound2 = np.log10(physiological_range(parnames[i])[1])

        params_run[parnames[i]] = 10**right_bound2
        sgs_run.set_parameters(params_run)
        if ((concept == Concept.OSCILLATION and sgs_run.is_oscillatory()) or (
                concept == Concept.BISTABILITY and sgs_run.is_bistable())):
            right_bound = right_bound2
        else:
            while right_bound2-right_bound1 > 0.01:
                mid = (right_bound2+right_bound1)/2
                params_run[parnames[i]] = 10**mid
                sgs_run.set_parameters(params_run)
                if ((concept == Concept.OSCILLATION and sgs_run.is_oscillatory()) or (
                        concept == Concept.BISTABILITY and sgs_run.is_bistable())):
                    right_bound1 = mid
                else:
                    right_bound2 = mid
            right_bound = (right_bound2+right_bound1)/2

        bounds[i,0] = np.log10(params[parnames[i]])
        bounds[i,1] = left_bound
        bounds[i,2] = right_bound

    # Define column header.
    columns = ['log10_var', 'log10_left_bound', 'log10_right_bound']

    # Save bifurcation data.
    df = pd.DataFrame(data = bounds, columns=columns, index=parnames)

    df.to_csv(fname)

def bifurcation_sampling(sgs, fname, bif_endpoints, concept = Concept.OSCILLATION):
    """ Check 5 points in bifurcation range for oscillations/bistability, in csv file.

    Bounds are limited by physiological range.
    """

    # -------------------------------------------------------------------------------
    # bifurcation parameters
    # -------------------------------------------------------------------------------

    # List all bifurcation parameters.
    if isinstance(sgs, DDDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3',
                   'f12', 'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3',
                   'bcoop12', 'bcoop23', 'bcoop13', 'bcoop123', 'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']
    elif isinstance(sgs, DDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12',
               'kb1', 'kb2', 'ku1', 'ku2', 'bcoop', 'ucoop']
    elif isinstance(sgs, MDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm',
               'kbd']

    # -------------------------------------------------------------------------------
    # bifurcation for every parameter
    # -------------------------------------------------------------------------------

    params = sgs.get_parameters()

    # Number of sets that will be checked per parameter
    N = 5

    bifs = np.zeros((len(parnames),N*3))

    sgs_run = type(sgs)()

    for i, par in enumerate(parnames):
        # Set parameters of "running" sgs to copy of original parameters.
        params_run = { x : params[x] for x in parnames }
        for x in ['taum', 'taumRNA', 'DNAtot']:
            params_run[x] = params[x]

        # -------------------------------------------------------------------------------
        # check points distributed in bifurcation range
        # -------------------------------------------------------------------------------

        # Bounds of the oscillating/bistable range.
        left_bound = bif_endpoints.loc[par,'log10_left_bound']
        right_bound = bif_endpoints.loc[par,'log10_right_bound']

        for j in range(N):
            params_run[par] = 10**(left_bound + (j+1)/(N+1)*(right_bound - left_bound))

            sgs_run.set_parameters(params_run)

            if concept == Concept.OSCILLATION:
                bif_run = sgs_run.oscillation_parameters('d')

                bifs[i,j*3+0] = bif_run[0]
                bifs[i,j*3+1] = bif_run[1]
                bifs[i,j*3+2] = bif_run[2]
            elif concept == Concept.BISTABILITY:
                bif_run = np.sort(sgs_run.steady_state()['d'])

                if len(bif_run) == 3:
                    bifs[i, j * 3 + 0] = bif_run[0]
                    bifs[i, j * 3 + 1] = bif_run[1]
                    bifs[i, j * 3 + 2] = bif_run[2]
                else:
                    for k in range(3):
                        bifs[i, j * 3 + k] = np.nan

    # Define column headers.
    columns = []
    for j in range(N):
        if concept == Concept.OSCILLATION:
            columns += ['p%d_per'%(j+1), 'p%d_min'%(j+1), 'p%d_max'%(j+1)]
        elif concept == Concept.BISTABILITY:
            columns += ['L%d_per'%(j+1), 'I%d_min'%(j+1), 'H%d_max'%(j+1)]

    # Save data.
    df = pd.DataFrame(data = bifs, columns=columns, index=parnames)

    df.to_csv(fname, float_format='%.3E')

def bifurcation_graph_data(sgs, fname, bif_endpoints, concept = Concept.OSCILLATION):
    """ Find and save oscillation/bistability parameters for 50 points in bifurcation range.

    Bounds are limited by physiological range.
    """

    # -------------------------------------------------------------------------------
    # bifurcation parameters
    # -------------------------------------------------------------------------------

    # List all bifurcation parameters.
    if isinstance(sgs, DDDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3',
                   'f12', 'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3',
                   'bcoop12', 'bcoop23', 'bcoop13', 'bcoop123', 'ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']
    elif isinstance(sgs, DDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12',
               'kb1', 'kb2', 'ku1', 'ku2', 'bcoop', 'ucoop']
    elif isinstance(sgs, MDS):
        parnames = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm',
               'kbd']

    # -------------------------------------------------------------------------------
    # bifurcation for every parameter
    # -------------------------------------------------------------------------------

    params = sgs.get_parameters()

    # Number of sets that will be checked per parameter
    N = 50

    bifs = np.zeros((len(parnames),N*4))

    sgs_run = type(sgs)()

    for i, par in enumerate(parnames):
        # Set parameters of "running" sgs to copy of original parameters.
        params_run = { x : params[x] for x in parnames }
        for x in ['taum', 'taumRNA', 'DNAtot']:
            params_run[x] = params[x]

        # -------------------------------------------------------------------------------
        # check points distributed in bifurcation range
        # -------------------------------------------------------------------------------

        # Bounds of the oscillating/bistable range.
        left_bound = bif_endpoints.loc[par,'log10_left_bound']
        right_bound = bif_endpoints.loc[par,'log10_right_bound']

        # Bounds physiological range.

        phys_left_bound = np.log10(physiological_range(par)[0])
        phys_right_bound = np.log10(physiological_range(par)[1])


        # Define points in and around bifurcation range to study.
        xl = np.zeros(0)
        xr = np.zeros(0)

        if left_bound != phys_left_bound:
            xl = np.append(xl, np.linspace(max(phys_left_bound, 0.7*left_bound-0.3*right_bound),
                            max(phys_left_bound, 0.9*left_bound-0.1*right_bound), 5, endpoint=False))
            xl = np.append(xl, np.linspace(max(phys_left_bound, 0.9 * left_bound - 0.1 * right_bound),
                                         max(phys_left_bound, 1.1 * left_bound + 0.1 * right_bound), 10, endpoint=False))
            xl = np.unique(xl)
        else:
            xl = np.append(xl, phys_left_bound)
        if right_bound != phys_right_bound:
            xr = np.append(xr, np.linspace(min(phys_right_bound, 0.9 * right_bound + 0.1 * left_bound),
                                         min(phys_right_bound, 1.1 * right_bound - 0.1 * left_bound), 11)[1:])
            xr = np.append(xr, np.linspace(min(phys_right_bound, 1.1 * right_bound - 0.1 * left_bound),
                                            min(phys_right_bound, 1.3 * right_bound - 0.3 * left_bound), 6)[1:])
            xr = np.unique(xr)
        else:
            xr = np.append(xr, phys_right_bound)

        xm = np.linspace(xl[-1], xr[0], N - len(xl) - len(xr) + 2)[1:-1]

        x = np.append(np.append(xl, xm), xr)

        for j, exp in enumerate(x):
            params_run[par] = 10**exp

            sgs_run.set_parameters(params_run)

            if concept == Concept.OSCILLATION:
                bif_run = sgs_run.oscillation_parameters('d')

                bifs[i,j*4+0] = exp
                bifs[i,j*4+1] = bif_run[0]
                bifs[i,j*4+2] = bif_run[1]
                bifs[i,j*4+3] = bif_run[2]
            elif concept == Concept.BISTABILITY:
                bif_run = np.sort(sgs_run.steady_state()['d'])

                if len(bif_run) == 3:
                    bifs[i, j * 4 + 0] = exp
                    bifs[i, j * 4 + 1] = bif_run[0]
                    bifs[i, j * 4 + 2] = bif_run[1]
                    bifs[i, j * 4 + 3] = bif_run[2]
                else:
                    for k in range(3):
                        bifs[i, j * 3 + k] = np.nan

    # Define column headers.
    columns = []
    for j in range(N):
        if concept == Concept.OSCILLATION:
            columns += ['p%d_log10(par)'%(j+1), 'p%d_per'%(j+1), 'p%d_min'%(j+1), 'p%d_max'%(j+1)]
        elif concept == Concept.BISTABILITY:
            columns += ['p%d_log10(par)'%(j+1), 'L%d_per'%(j+1), 'I%d_min'%(j+1), 'H%d_max'%(j+1)]

    # Save data.
    df = pd.DataFrame(data = bifs, columns=columns, index=parnames)

    df.to_csv(fname, float_format='%.3E')

def bifurcation_analysis(sgs, system, concept):
    """ Find bifurcation points around sgs solution and do small sampling for continuity."""

    folder = foldername(sgs, system, concept)

    df = pd.read_csv(folder + "variables.csv", index_col=0)

    if concept == Concept.OSCILLATION:
        list_indices = df[~np.isnan(df['period'])].index.tolist()
    elif concept == Concept.BISTABILITY:
        list_indices = df.index.tolist()

    one_core_ = partial(one_core, sgs=sgs, folder=folder, concept=concept)

    # one_core_(list_indices[0])

    p = multiprocessing.Pool(num_cores)
    p.map(one_core_, list_indices)

def oscillation_discontinuity(sgs, system, concept):
    """ Read all bifurcation files and print how many are discontinuous."""

    folder = foldername(sgs, system, concept)

    df = pd.read_csv(folder + "variables.csv", index_col=0)

    list_indices = df[~np.isnan(df['period'])].index.tolist()

    count_nan = 0
    count_stab = 0

    for idx in list_indices:
        bifpoints_name = folder + "bifurcations/points-%d.csv" % idx
        bif_name = folder + "bifurcations/bifurcation-%d.csv" % idx

        if pyos.path.exists(bif_name):
            data = pd.read_csv(bif_name, index_col=0)

            if data.isnull().values.any():
                print(bif_name, sum(len(data) - data.count())/3)

                count_nan += 1
            else:
                count_stab += 1

    print("nan, stab", count_nan, count_stab, count_nan/count_stab)

    return

def main():
    sgs = MDS
    system = System.RANDOM
    concept = Concept.BISTABILITY

    bifurcation_analysis(sgs, system, concept)

    #oscillation_discontinuity(sgs, system, concept)

    #filename = "oscillation/2DS/variables.csv"
    #data = pd.read_csv(filename, index_col=0)

    #params = data.loc[2].to_dict()
    #os = DDS(params)

    #bif_endpoints = pd.read_csv("oscillation/2DS/bifurcations/points-2.csv", index_col=0)
    #bifurcation_graph_data(os, "oscillation/2DS/bifurcations/bif-complete2.csv", bif_endpoints, Concept.OSCILLATION)

if __name__ == '__main__':
    main()