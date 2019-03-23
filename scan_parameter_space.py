#!/usr/bin/env python

"""scan.py: Scan parameter space for oscillations/bistability (O/B).
            Random scanning for O and B, and systematic scanning for B."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from SGS import *
import multiprocessing
from functools import partial
import time
import itertools as it

oscillator = DDDS

num_cores = 2


# --------------------------------------------------------------------------------------------------------
# RANDOM SCAN
# --------------------------------------------------------------------------------------------------------

def one_random(r, fpar, sgs, system=System.RANDOM, concept=Concept.OSCILLATION):
    """ Check one random set of parameters on oscillatory behavior/bistability.
    If oscillations/bistability, save set in fpar.

    r is random seed """

    # Generate a system with random parameters.
    found_good_parameters = False
    tries = 0
    while (not found_good_parameters):
        pars = random_params(sgs, physiological_range, system=system, seed=int(time.time() % r * 1.5e9 + r * tries))

        s = sgs()
        s.set_parameters(pars);

        if system != System.SSLRPB or s.SsLrpB_compatible():
            found_good_parameters = True

        tries += 1

    # Check if oscillatory/bistable and save if so.
    if concept == Concept.OSCILLATION and s.isoscillatory():
        with open(fpar, "a") as file:
            file.write(s.parameter_line(sep=","))
    elif concept == Concept.BISTABILITY:
        ss = np.sort(s.steady_state()['d'])

        if (s.is_bistable() and ss[1] - ss[0] > 10 and 500 < ss[2] < 600 and ss[1] < 100):
            with open(fpar, "a") as file:
                file.write(s.parameter_line(sep=","))

    # Delete the SGS.
    del s


def random_search(fpar, sgs, system=System.RANDOM, concept=Concept.OSCILLATION):
    """ Check 2e6 random parameter sets for oscillations. """

    # Write names in parameters file
    l = ''
    for p in sgs.PARAMETERS:
        l += '%s,' % p
    l = l[:-1] + '\n'

    with open(fpar, "a") as file:
        file.write(l)

    # Set multiprocessing.
    p = multiprocessing.Pool(num_cores)

    np.random.seed(int(time.time()))
    onerandom_ = partial(one_random, fpar=fpar, sgs=sgs, system=system, concept=concept)
    for i in range(10):
        # Set random number before start of multiprocessing to avoid that different cores use the same random sets.
        r = [np.random.uniform(0, 1) for _ in range(int(2e5))]
        p.map(onerandom_, r)
    return

# --------------------------------------------------------------------------------------------------------
# OSCILLATIONS
# --------------------------------------------------------------------------------------------------------

def add_oscillation_data(file, oss):
    """ Add parameters with respect to oscillation (period, minima and maxima of variables) to dataframe.

    oss is the type of oscillator (MDS, DDS or DDDS). """

    # Set variables of which minima and maxima will be determined.
    if oss == MDS:
        vars = ['DNA0', 'DNAm', 'DNAd', 'mRNA', 'm', 'd']
    elif oss == DDS:
        vars = ['DNA0', 'DNA1', 'DNA2', 'DNA12', 'mRNA', 'm', 'd']
    elif oss == DDDS:
        vars = ['DNA0', 'DNA1', 'DNA2', 'DNA3', 'DNA12', 'DNA13', 'DNA23', 'DNA123', 'mRNA', 'm', 'd']

    # Read data.
    df = pd.read_csv(file)

    # Add oscillation parameters in columns to dataframe.
    df['period'] = np.nan
    for var in oss.ALL_VARS:
        df['%s_min' % var] = np.nan;
        df['%s_max' % var] = np.nan

    os = oss()
    for index, row in df.iterrows():
        params = row.to_dict()  # readparams(i, data, oss)
        os.set_parameters(params)
        oscillation_parameters = os.oscillation_parameters(oss.ALL_VARS)

        # Add period, minima and maxima to dataframe.
        df.loc[index, 'period'] = oscillation_parameters[0]
        for i, var in enumerate(oss.ALL_VARS):
            df.loc[index, '%s_min' % var] = oscillation_parameters[1 + 2 * i]
            df.loc[index, '%s_max' % var] = oscillation_parameters[2 + 2 * i]

    # Save dataframe.
    df.to_csv("oscillation/MDS/test.csv")


def meets_selection_criteria(mRNAmax, mmax, dmax, mRNAamp, mamp, damp):
    """ Return whether the variable meet the selection criteria. """

    if mRNAmax < 20 and mmax < 5000 and dmax < 5000 and mRNAamp > 1.0 and mamp > 1.0 and damp > 10.0 and mRNAamp / mRNAmax > 0.33 and mamp / mmax > 0.33 and damp / dmax > 0.33:  # and dmax > 2*mmax:
        return True
    return False


def add_selection_monotonicity_data(file, oss):
    """ Add monotonicity type and whether solution passes selection criteria to dataframe of file.

    oss is the type of oscillator (MDS, DDS or DDDS). """

    # Read data.
    df = pd.read_csv(file, index_col=0)

    os = oss()

    # Add oscillation parameters in columns to dataframe.
    df['selected'] = np.nan
    df['monotonicity_type'] = np.nan

    for index, row in df.iterrows():
        params = row.to_dict()
        os.set_parameters(params)

        # Add period, minima and maxima to dataframe.
        df.loc[index, 'monotonicity_type'] = os.monotonicity_type()
        df.loc[index, 'selected'] = meets_selection_criteria(row['mRNA_max'], row['m_max'], row['d_max'],
            (row['mRNA_max'] - row['mRNA_min']), (row['m_max'] - row['m_min']), (row['d_max'] - row['d_min']))

    # Save file.
    df.to_csv(file)

# --------------------------------------------------------------------------------------------------------
# BISTABILITY SYSTEMATIC SCAN
# --------------------------------------------------------------------------------------------------------

def check_SsLrpB_compatibility():
    data = pd.read_csv("bistability/SsLrpB/variables.csv", index_col=0)

    list_indices = data.index.tolist()

    sgs = DDDS()

    print("Solution is SsLprB compatible (True) : mixed feedback")

    for i in list_indices:
        params = data.loc[i].to_dict()

        sgs.set_parameters(params)

        print(i, sgs.SsLrpB_compatible())

def induction_time(sgs, p, Q):
    ''' Return induction time from the intermediate steady state to high steady state of the bistable switch.'''

    sg = sgs()

    # determine parameters from A-F parameters
    if sgs == DDS:
        a, b, c, d = p

        f12 = a / c

        Kd1 = np.random.uniform(0, 1)
        Kd2 = np.random.uniform(0, 1)
        s = Kd1 + Kd2
        Kd1 *= d / s
        Kd2 *= d / s
        omega = c / (Kd1 * Kd2)
        bcoop = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega)), min(2, 2 + np.log10(omega)))
        ucoop = bcoop / omega
        f1 = 10 ** np.random.uniform(max(-3, np.log10((b - Kd2 * 1e2) / Kd1)), min(2, np.log10((b - Kd2 * 1e-3) / Kd1)))
        f2 = (b - f1 * Kd1) / Kd2
        beta = 0.1 * 0.01 * np.sqrt((2 + 0.001) / 0.01) / (
            Q * 1)  # gammam*gammamRNA*sqrt((alphadiss+gammad)/alphaass)/(Q*phi0)

        params = {'beta': beta, 'gammam': 0.1, 'gammamRNA': 0.01, 'gammad': 0.001, 'alphaass': 0.01, 'alphadiss': 2,
                  'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'phi0': 1, 'f1': f1, 'f2': f2, 'f12': f12, 'kb1': 1, 'Kd1': Kd1,
                  'kb2': 1, 'Kd2': Kd2, 'bcoop': bcoop, 'ucoop': ucoop, 'vol': 1}

        sg.set_parameters(params)

    elif sgs == MDS:
        Km, Kd, fm, fd = p

        gammam = 0.1
        gammamRNA = 0.01
        phi0 = 1
        beta = gammam * gammamRNA / (Q * phi0)
        alphaass = 0.01
        alphadiss = 2
        gammad = 0.001
        params = {'beta': beta, 'gammam': 0.1, 'gammamRNA': 0.01, 'gammad': gammad, 'alphaass': alphaass,
                  'alphadiss': alphadiss, 'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'phi0': 1, 'fm': fm, 'fd': fd,
                  'kbm': 1, 'Km': Km, 'kbd': 1, 'Kd': Kd * (alphadiss + gammad) / alphaass, 'vol': 1}

    elif sgs == DDDS:
        a, b, c, d, e, f = p

        f123 = a / d

        x = np.linspace(0, f, 10000)

        found = False
        if f < 3.5e-4 or (1e-8 + 2e-4 * (f - 2e-4)) * 1e-4 > 0.9 * e:
            Kd1 = 1e-4
            Kd2 = 1e-4
            Kd3 = f - Kd1 - Kd2
            found = True
        elif f > 2.5e4:
            Kd1 = 1e4
            Kd2 = 1e4
            Kd3 = f - Kd1 - Kd2
            found = True
        elif f ** 2 / 3 * 1e3 < 1.1 * e:
            Kd1 = f / 3
            Kd2 = f / 3
            Kd3 = f - Kd1 - Kd2
            found = True
        while not found:
            #fKd3 = lambda x: -x ** 3 + f * x ** 2 - e / 1e3 * x + d / 1e10
            #zeroes = np.where(np.diff(np.sign(fKd3(x))))[0]
            #Kd3 = brentq(fKd3, 0, x[zeroes[0] + 1], xtol=1e-10, rtol=1e-5, maxiter=200)

            lb = max(1e-4, (f - 2 * np.sqrt(f ** 2 - 3 * e / 1e3)) / 3)
            rb = np.nanmin([f, (f - np.sqrt(f ** 2 - 4 * e / 1e3)) / 2])
            Kd1 = 10 ** np.random.uniform(np.log10(lb), np.log10(rb))
            Kd2 = (-(Kd1 - f) + np.sqrt((Kd1 - f) ** 2 - 4 * (1e-3 * e + Kd1 * (Kd1 - f)))) / 2
            Kd3 = (-(Kd1 - f) - np.sqrt((Kd1 - f) ** 2 - 4 * (1e-3 * e + Kd1 * (Kd1 - f)))) / 2
            if 1.1e4 >= Kd2 >= 0.9e-4 and 1.1e4 >= Kd3 >= 0.9e-4:
                found = True
            else:
                Kd1 = (f - np.sqrt(f ** 2 - 3 * e * 1e-3)) / 3
                Kd2 = Kd1
                Kd3 = f - Kd1 - Kd2
                found = True

        found = False
        mine = -3;
        maxe = 2
        while not found:
            f1 = 10 ** np.random.uniform(mine, maxe)
            f2 = 10 ** np.random.uniform(mine, maxe)
            f3 = 10 ** np.random.uniform(mine, maxe)
            C = f1 * Kd1 + f2 * Kd2 + f3 * Kd3
            f1 *= c / C;
            f2 *= c / C;
            f3 *= c / C
            if min([f1, f2, f3]) > 0.0009 and max([f1, f2, f3]) < 101:
                found = True
            else:
                mine += 0.5;
                maxe -= 0.5
                if mine > maxe:
                    mine = 0;
                    maxe = 0
                    found = True

        found = False
        mine = np.log10(0.05)
        maxe = np.log10(20)
        while not found:
            omega12 = 10 ** np.random.uniform(mine, maxe)
            omega13 = 10 ** np.random.uniform(mine, maxe)
            omega23 = 10 ** np.random.uniform(mine, maxe)
            E = omega12 * Kd1 * Kd2 + omega23 * Kd2 * Kd3 + omega13 * Kd1 * Kd3
            omega12 *= e / E
            omega23 *= e / E
            omega13 *= e / E
            if min([omega12, omega23, omega13]) > 0.049 and max([omega12, omega23, omega13]) < 20.1:
                found = True
            else:
                mine += 0.1
                maxe -= 0.1
                if mine > maxe:
                    mine = 0;
                    maxe = 0
                    found = True

        found = False
        mine = -3
        maxe = 2
        while not found:
            f12 = 10 ** np.random.uniform(mine, maxe)
            f23 = 10 ** np.random.uniform(mine, maxe)
            f13 = 10 ** np.random.uniform(mine, maxe)
            B = f12 * omega12 * Kd1 * Kd2 + f23 * omega23 * Kd2 * Kd3 + f13 * omega13 * Kd1 * Kd3
            f12 *= b / B
            f23 *= b / B
            f13 *= b / B
            if min([f1, f2, f3]) > 0.0009 and max([f1, f2, f3]) < 101:
                found = True
            else:
                mine += 0.5
                maxe -= 0.5
                if mine > maxe:
                    mine = 0
                    maxe = 0
                    found = True

        omega123 = d / (Kd1 * Kd2 * Kd3 * omega12 * omega23 * omega13)
        if omega123 > 20 or omega123 > 0.05:
            print('PROBLEM', omega123)

        bcoop12 = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega12)), min(2, 2 + np.log10(omega12)))
        ucoop12 = bcoop12 / omega12

        bcoop23 = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega23)), min(2, 2 + np.log10(omega23)))
        ucoop23 = bcoop23 / omega23

        bcoop13 = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega13)), min(2, 2 + np.log10(omega13)))
        ucoop13 = bcoop13 / omega13

        bcoop123 = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega123)), min(2, 2 + np.log10(omega123)))
        ucoop123 = bcoop123 / omega123

        beta = 0.1 * 0.01 * np.sqrt((1 + 0.001) / 0.01) / (Q * 1)
        # gammam*gammamRNA*sqrt((alphadiss+gammad)/alphaass)/(Q*phi0)

        params = {'beta': beta, 'gammam': 0.1, 'gammamRNA': 0.01, 'gammad': 0.001, 'alphaass': 0.01, 'alphadiss': 1,
                  'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'phi0': 1, 'f1': f1, 'f2': f2, 'f3': f3, 'f12': f12,
                  'f23': f23, 'f13': f13, 'f123': f123, 'kb1': 1, 'Kd1': Kd1, 'kb2': 1, 'Kd2': Kd2, 'kb3': 1,
                  'Kd3': Kd3, 'bcoop12': bcoop12, 'ucoop12': ucoop12, 'bcoop13': bcoop13, 'ucoop13': ucoop13,
                  'bcoop23': bcoop23, 'ucoop23': ucoop23, 'bcoop123': bcoop123, 'ucoop123': ucoop123, 'vol': 1}

    sg.set_parameters(params)

    return sg.deterministic_induction_time()

def monotonicity_type(sgs, p):
    ''' Return monotonicity type for A-F parameters.'''

    xx = np.linspace(0, 1000, 1000)

    # reconstruct response curve
    if sgs == DDS:
        a, b, c, d = p

        y = (a * xx + b * np.sqrt(xx) + 1) / (c * xx + d * np.sqrt(xx) + 1)
    elif sgs == MDS:
        Km, Kd, fm, fd = p

        y = (fd * Kd * xx ** 2 + fm * Km * xx + 1) / (Kd * xx ** 2 + Km * xx + 1)
    elif sgs == DDDS:
        a, b, c, d, e, f = p

        y = (a * xx ** 3 + b * xx ** 2 + c * xx + 1) / (d * xx ** 3 + e * xx ** 2 + f * xx + 1)

    # derivative of response curve
    dy = y[1:] - y[:-1]

    # check for local extrema (derivative changes sign)
    zero_crossings = np.where(np.diff(np.sign(dy)))[0]

    # determine non-monotonicity type
    if len(zero_crossings) == 1: # and zero_crossings[0] != 0):
        if dy[zero_crossings[0]] > 0: #dy[0] > 0 and dy[-1] < 0:
            nonmono = 1
        elif dy[zero_crossings[0]] < 0: #dy[0] < 0 and dy[-1] > 0:
            nonmono = 2
    elif len(zero_crossings) > 1:
        nonmono = 3
    else:
        if y[0]<y[-1]:
            nonmono = 0
        if y[0] > y[-1]:
            nonmono = 4
    if False: # check validity of function
        fig = plt.figure('%d'%nonmono)
        ax = fig.add_subplot(211)
        ax.plot(xx, y)
        ax = fig.add_subplot(212)
        ax.plot(xx[:-1], dy)
        plt.show()

    return nonmono

def check_bistability(p=[], file=0, Hs=[], sgs=MDS):
    ''' Check if response function dictated by parameters p can lead to bistability
    with high steady state in interval Hs, save in file if it is the case '''

    # protein concentration (monomer for MDS, dimer for 2DS and 3DS)
    x = np.linspace(0, 1000, 10000)

    if sgs == MDS:
        Km, Kd, fm, fd = p
        dx = np.append(x[1:] - x[:-1], 0)

        def response(x):
            return (fd * Kd * x ** 2 + fm * Km * x + 1) / (Kd * x ** 2 + Km * x + 1)

        def derivative_response(x):
            return (2 * x * (fd * Kd - Kd) + Km * (fd * Kd * x ** 2 - 1) - fm * Km * Kd * x ** 2 + fm * Km) / (
                    Kd * x ** 2 + Km * x + 1) ** 2

        def dxdt(x, q):
            return response(x) - q * x

        Q = np.linspace(response(Hs[0]) / Hs[0], response(Hs[1]) / Hs[1], 50);
    elif sgs == DDS:
        a, b, c, d = p
        dx = x * np.append(x[1:] - x[:-1], 0) # monomer dx

        def response(x):
            return (a * x ** 2 + b * x + 1) / (c * x ** 2 + d * x + 1)

        def derivative_response(x):
            return (2 * x * (a - c) + d * (a * x ** 2 - 1) - b * c * x ** 2 + b) / (Kd * x ** 2 + Km * x + 1) ** 2

        def dxdt(x, q):
            return response(x) - q * np.sqrt(x)

        Q = np.linspace(response(Hs[0]) / np.sqrt(Hs[0]), response(Hs[1]) / np.sqrt(Hs[1]), 50)
    elif sgs == DDDS:
        a, b, c, d, e, f = p

        x = np.linspace(0, 1000, 10000)  # dimer spacing
        dx = x * np.append(x[1:] - x[:-1], 0)  # monomer dx

        def response(x):
            return (a * x ** 3 + b * x ** 2 + c * x + 1) / (d * x ** 3 + e * x ** 2 + f * x + 1)

        def derivative_response(x):
            return ((3 * a * x ** 2 + 2 * b * x + c) * (d * x ** 3 + e * x ** 2 + f * x + 1) - (
                    3 * d * x ** 2 + 2 * e * x + f) * (a * x ** 3 + b * x ** 2 + c * x + 1)) / (
                    d * x ** 3 + e * x ** 2 + f * x + 1) ** 2

        def dxdt(x, q):
            return response(x) - q * np.sqrt(x)


        Q = np.linspace(response(Hs[0]) / np.sqrt(Hs[0]), response(Hs[1]) / np.sqrt(Hs[1]), 50)
    s = [np.nan] * 100

    for m, q in enumerate(Q):
        y = dxdt(x, q)
        zero_crossings = np.where(np.diff(np.sign(y)))[0]
        Nss = len(zero_crossings)
        if (Nss > 2):
            L = x[zero_crossings[0]]
            I = x[zero_crossings[1]]
            if (I < 100 and I - L > 10):
                s[m] = sum((dx / y)[np.logical_and(y > 0, x > I)])

    if (not np.all(np.isnan(s))):
        # Find value with optimal score
        Q = Q[np.nanargmin(s)]

        # Determine monotonicity type
        nonmono = monotonicity_type(sgs, p)

        # Determine steady states and score
        zero_crossings = np.where(np.diff(np.sign(dxdt(x, Q))))[0]
        L = x[zero_crossings[0]]
        I = x[zero_crossings[1]]
        H = x[zero_crossings[2]]
        R = sum((dx / y)[np.logical_and(y > 0, x > I)])

        # Calculate induction time
        T = induction_time(sgs, p, Q)

        # Save values
        if file != 0:
            s = ""
            for x in [R, H, I, L, Q]:
                s += "%.5E," % x
            for x in p:
                s += "%.5E," % x
            s += "%d," % nonmono
            s += "%.5E\n" % T
            with open(file, 'a') as f:
                f.write(s)
        return [R, H, I, L, Q, nonmono, T]
    else:
        return []

def bistability_scan(sgs, Hs):
    if sgs == DDS:
        file = 'results/2DSH%d-%d-newdef.txt' % (Hs[0], Hs[1])

        NN = 20

        Ds = np.logspace(np.log10(2e-4), np.log10(2e4), NN)  # Kd1 + Kd2

        with open(file, 'a') as f:
            f.write('R,H,I,L,Q,A,B,C,D,nonmono,T\n')

        for i, d in enumerate(Ds):
            Bs = np.logspace(-3, 2, NN) * d
            for j, b in enumerate(Bs):
                print(i, j)
                Cs = np.logspace(-4 + np.log10(d - 1e-4) - 4, 2 * np.log10(d / 2) + 4, NN)
                for c in Cs:
                    As = np.logspace(-3 + np.log10(c), 2 + np.log10(c), NN)
                    for a in As:
                        check_bistability([a, b, c, d], file, Hs, DDS)

    elif sgs == MDS:
        file = 'results/MDSH%d-%d.txt' % (Hs[0], Hs[1])
        NN = 30
        Kms = np.logspace(-7, 3, NN);  # -7 -> 3
        Kds = np.logspace(-10, -3, NN);  # -10 -> 3 Kd*alphaass/(alphadiss + gammad)
        fms = np.logspace(-3, 2, NN);  # -3 -> 2
        fds = np.logspace(-3, 2, NN)  # -3 -> 2

        fs = tuple(it.product(fms, fds))

        p = multiprocessing.Pool(num_cores)

        check_bistability_MDS = partial(check_bistability, file=file, Hs=Hs, sgs=MDS)

        with open(file, 'a') as f:
            f.write('R,H,I,L,Q,Km,Kd,fm,fd,nonmono,T\n')
        for Km in Kms:
            for Kd in Kds:
                set = [[Kd, Km, f[0], f[1]] for f in fs]
                p.map(check_bistability_MDS, set)

    elif sgs == DDDS:
        file = 'bistability/3DS/test.csv'  # H%d-%d.csv' % (Hs[0], Hs[1])

        NN = 8

        with open(file, 'a') as f:
            f.write('R,H,I,L,Q,A,B,C,D,E,F,nonmono,T\n')

        p = multiprocessing.Pool(num_cores)

        # Loop over all values of parameters A-E and look for bistable solutions
        for f in np.logspace(np.log10(3e-4), np.log10(3e4), NN):
            # To find extreme values of the parameter, load hull information
            data = np.loadtxt('results/hull/%.3E.txt' % f)
            if (f - 2e-4 <= 1e4):
                minx = 1e-8 + 2 * 1e-4 * (f - 2e-4)
            elif (f - 1e-4 - 1e4 <= 1e4):
                minx = (1e-4 * (f - 1e4 - 1e-4) + 1e4 * (f - 1e4 - 1e-4) + 1)
            else:
                minx = ((f - 2e4) ** 2 + 1e4 * (f - 2e4) + 1e8)
            for e in np.logspace(np.log10(0.05 / 20) + np.log10(minx), np.log10(f ** 2 / 3) + np.log10(20 / 0.05), NN):
                minD = 0;
                maxD = 0;
                for i in range(len(data)):
                    if (0.9 * e < data[i, 0] < 1.1 * e):
                        _, minD, minKd1, minKd2, minKd3, minomega12, minomega23, minomega13, maxD, maxKd1, maxKd2, maxKd3, maxomega12, maxomega23, maxomega13 = \
                            data[i]
                for d in np.logspace(np.log10(minD), np.log10(maxD), NN):
                    As = np.logspace(-3 + np.log10(d), 2 + np.log10(d), NN)
                    Bs = np.logspace(-3 + np.log10(e), 2 + np.log10(e), NN)
                    Cs = np.logspace(-3 + np.log10(f), 2 + np.log10(f), NN)

                    check_bistability_ = partial(check_bistability, file=file, Hs=Hs, sgs=DDDS)
                    values = tuple([[a, b, c, d, e, f] for a in As for b in Bs for c in Cs])
                    p.map(check_bistability_, values)

def check_bistability_SsLrpB(f13=0, f123=0):
    ''' Save the score, high steady state, degradation rate, C, B and induction time of SsLrpB systems
    with optimized induction time'''
    print(f13, f123)

    # Make 3DS
    ddds = DDDS()

    # Set fixed parameters for SsLrpB compatible system
    ct = 1e-3 / 2.4

    Kd1 = 73.5 * ct
    Kd2 = 0.7 * ct
    Kd3 = 49.1 * ct
    omega12 = 4.1
    omega13 = 7.9
    omega23 = 2.1
    omega123 = 3.1

    # Calculate parameters D, E and F
    d = omega12 * omega23 * omega13 * omega123 * Kd1 * Kd2 * Kd3
    e = omega12 * Kd1 * Kd2 + omega23 * Kd2 * Kd3 + omega13 * Kd1 * Kd3
    f = Kd1 + Kd2 + Kd3

    # Set ranges for A, B and C parameters
    Cs = np.logspace(-3 + np.log10(Kd1 + Kd2 + Kd3), 2 + np.log10(Kd1 + Kd2 + Kd3), 50)  # 50);
    a = 10 ** f123 * omega12 * omega23 * omega13 * omega123 * Kd1 * Kd2 * Kd3;
    Bs = 10 ** f13 * omega13 * Kd1 * Kd3 + 10 ** np.linspace(-3, 2, 20) * (omega12 * Kd1 * Kd2 + omega23 * Kd2 * Kd3);

    # Generate name to save file (f123, f13 * 10)
    sf123 = int(round(f123 * 10));
    if (sf123 < 0):
        sf123 = 'm%d' % -sf123
    else:
        sf123 = '%d' % sf123
    sf13 = int(round(f13 * 10));
    if (sf13 < 0):
        sf13 = 'm%d' % -sf13
    else:
        sf13 = '%d' % sf13

    # Allocate memory to store score (R), induction time (T), high steady state (H) and 'degragadation rate' (Q)
    R = np.empty([len(Cs), len(Bs)]);
    R[:] = np.nan;
    T = np.empty([len(Cs), len(Bs)]);
    T[:] = np.nan;
    H = np.empty([len(Cs), len(Bs)]);
    H[:] = np.nan;
    Q = np.empty([len(Cs), len(Bs)]);
    Q[:] = np.nan;

    # Number of Qs to check for bistability
    N = 100

    # Loop through all B and C values
    Bs, Cs = np.meshgrid(Bs, Cs)
    for i in range(len(Bs)):
        for j in range(len(Bs[0])):
            c = Cs[i, j];
            b = Bs[i, j];

            x = np.linspace(0, 1000, 10000);
            dx = x * np.append(x[1:] - x[:-1], 0)  # monomer dx

            fy = lambda t: (a * t ** 3 + b * t ** 2 + c * t + 1) / (d * t ** 3 + e * t ** 2 + f * t + 1);

            # Set Q values
            sqrtc = np.linspace(fy(200) / np.sqrt(200), fy(1000) / np.sqrt(1000), N);

            # Allocate space
            crosl = [0 for _ in range(N)];
            crosm = [0 for _ in range(N)];
            crosr = [0 for _ in range(N)];
            s = [0 for _ in range(N)];

            for l in range(len(sqrtc)):
                # Calculate response function given q
                y = fy(x) - sqrtc[l] * np.sqrt(x);

                # Check if response is bistable (more than 2 zeros)
                zero_crossings = np.where(np.diff(np.sign(y)))[0]
                Nss = len(zero_crossings)
                if (Nss > 2):
                    crosl[l] = zero_crossings[0]
                    crosm[l] = zero_crossings[1]
                    crosr[l] = zero_crossings[2]
                    # If low and intermediate steady state are far enough apart and intermediate state is not too high
                    # (I - L > 10 and I < 100), calculate the score (integral of dx/dt between I and H)
                    if ((x[crosm[l]] - x[crosl[l]]) > 10 and x[crosm[l]] < 100):
                        s[l] = sum((dx / y)[np.logical_and(y > 0, x > x[crosm[l]])])

            # If any good bistable system was found (any score was given)
            if (max(s) > 0):
                # Determine best solution
                idxmax = s.index(max(s))

                # Set new Q ranges around best solution
                sqrtc2 = np.linspace(sqrtc[max(0, idxmax - 1)], sqrtc[min(idxmax + 2, N)], N);

                x2 = np.linspace(0, x[min(3 * crosr[idxmax], len(x) - 1)], 10000)
                dx2 = x2 * np.append(x2[1:] - x2[:-1], 0)  # monomer dx

                # Calculate new scores to find best one
                s = [0 for _ in range(N)]
                for l in range(N):
                    y = fy(x2) - sqrtc2[l] * np.sqrt(x2);
                    zero_crossings = np.where(np.diff(np.sign(y)))[0]
                    Nss = len(zero_crossings)
                    if (Nss > 2):
                        crosl[l] = zero_crossings[0]
                        crosm[l] = zero_crossings[1]
                        crosr[l] = zero_crossings[2]
                        if ((x2[crosm[l]] - x2[crosl[l]]) > 10 and x2[crosm[l]] < 100):
                            s[l] = sum((dx2 / y)[np.logical_and(y > 0, x2 > x2[crosm[l]])])

                # Save best values
                idxmax = s.index(max(s));
                H[i, j] = x2[crosr[idxmax]];
                R[i, j] = s[idxmax];
                Q[i, j] = sqrtc2[idxmax];

                # Generate random f1, f2 and f3 values and bcoop and ucoop values that result in the given A-F parameters
                found = False
                mine = -3;
                maxe = 2
                while (not found):
                    f1 = 10 ** np.random.uniform(mine, maxe)
                    f2 = 10 ** np.random.uniform(mine, maxe)
                    f3 = 10 ** np.random.uniform(mine, maxe)
                    C = f1 * Kd1 + f2 * Kd2 + f3 * Kd3
                    f1 *= c / C;
                    f2 *= c / C;
                    f3 *= c / C;
                    if (min([f1, f2, f3]) > 0.0009 and max([f1, f2, f3]) < 101):
                        found = True
                    else:
                        mine += 0.5;
                        maxe -= 0.5;
                        if (mine > maxe):
                            mine = 0;
                            maxe = 0;
                            found = True;

                bcoop12 = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega12)), min(2, 2 + np.log10(omega12)))
                ucoop12 = bcoop12 / omega12

                bcoop23 = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega23)), min(2, 2 + np.log10(omega23)))
                ucoop23 = bcoop23 / omega23

                bcoop13 = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega13)), min(2, 2 + np.log10(omega13)))
                ucoop13 = bcoop13 / omega13

                bcoop123 = 10 ** np.random.uniform(max(-3, -2 + np.log10(omega123)), min(2, 2 + np.log10(omega123)))
                ucoop123 = bcoop123 / omega123

                beta = 0.1 * 0.01 * np.sqrt((2 + 0.001) / 0.01) / (Q[i, j] * 1)

                f12 = (b - 10 ** f13 * omega13 * Kd1 * Kd3) / (omega12 * Kd1 * Kd2 + omega23 * Kd2 * Kd3);
                f23 = f12

                # Fix other parameters

                params = {'beta': beta, 'gammam': 0.1, 'gammamRNA': 0.01, 'gammad': 0.001, 'alphaass': 0.01,
                          'alphadiss': 2, 'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'phi0': 1, 'f1': f1, 'f2': f2, 'f3': f3,
                          'f12': f12, 'f23': f23, 'f13': 10 ** f13, 'f123': 10 ** f123, 'kb1': 1, 'Kd1': Kd1, 'kb2': 1,
                          'Kd2': Kd2, 'kb3': 1, 'Kd3': Kd3, 'bcoop12': bcoop12, 'ucoop12': ucoop12, 'bcoop13': bcoop13,
                          'ucoop13': ucoop13, 'bcoop23': bcoop23, 'ucoop23': ucoop23, 'bcoop123': bcoop123,
                          'ucoop123': ucoop123, 'vol': 1}

                # Do deterministic timeseries and determine the induction time
                ddds.setparameters(params)
                T[i, j] = ddds.parambistab();
            else:
                H[i, j] = np.nan;
                R[i, j] = np.nan;

        # Save all values
        data = np.vstack((R.flatten(), H.flatten(), Q.flatten(), Cs.flatten(), Bs.flatten(), T.flatten()))

        np.savetxt('new/f123-%s-f13-%s.txt' % (sf123, sf13), data.T, fmt='%.3E')

def bistability_scan_SsLrpB():
    ''' Save all information about bistability in the parameter space of the SsLrpB system'''

    p = multiprocessing.Pool(num_cores)

    # scan over all values of f123 and f13
    f13s = np.linspace(-3, 2, 51);
    f123s = np.linspace(-3, 2, 51);
    f123s = f123s[f123s >= -0.2]
    for f123 in f123s:
        check_bistability_SsLrpB_ = partial(check_bistability_SsLrpB, f123=f123)
        p.map(check_bistability_SsLrpB_, f13s)

# --------------------------------------------------------------------------------------------------------
# COUNT SOLUTIONS
# --------------------------------------------------------------------------------------------------------

def count_solutions_oscillations(file):
    """ Count the different types of monotonicity and print the results."""

    data = pd.read_csv(file)

    print("Without selection criteria",
          len(data[np.logical_and(np.logical_not(np.isnan(data['period'])), data['monotonicity_type'] == 0)]),
          len(data[np.logical_and(np.logical_not(np.isnan(data['period'])), data['monotonicity_type'] == 1)]),
          len(data[np.logical_and(np.logical_not(np.isnan(data['period'])), data['monotonicity_type'] == 2)]),
          len(data[np.logical_and(np.logical_not(np.isnan(data['period'])), data['monotonicity_type'] == 3)]),
          len(data[~np.isnan(data['period'])]))

    print("With selection criteria", len(data[np.logical_and(data['selected'], data['monotonicity_type'] == 0)]),
          len(data[np.logical_and(data['selected'], data['monotonicity_type'] == 1)]),
          len(data[np.logical_and(data['selected'], data['monotonicity_type'] == 2)]),
          len(data[np.logical_and(data['selected'], data['monotonicity_type'] == 3)]), len(data[data['selected']]))

def count_solutions_bistability():
    systems = ['MDS', '2DS', '3DS']

    for system in systems:
        f = 'bistability/%s/H500-600.csv' % system
        data = pd.read_csv(f, na_values='NAN')
        data = data[np.logical_and(data['nonmono'] != 4, np.log10(data['T']) > 0)]
        data.dropna(inplace=True)

        print(system, data.groupby('nonmono').count()['R'])
        print("total %s" % system, len(data))


def main():
    #for jj in range(2):
    #    random_search('test-MDS.csv', sgs=MDS, system=System.RANDOM, concept=Concept.BISTABILITY)
    # check_SsLrpB_compatibility()
    return

if __name__ == '__main__':
    main()

