"""scan_random_stoch.py"""

__author__ = "Sophie de Buyl"
__email__ = "Sophie.de.Buyl@vub.be"

# motivation: it seems that the vast majority of parameter sets leads the system to reach steady-state. we wonder
# (a) if a deterministic steady-state gives rise to a bursting behavior stochastically
# (b) if the noise is increasing/decreasing with more boxes added. 
# it seems more likely that the complexity will reduce the noise. adding boxes may help finetuning in the steady state level
# idee: look at noise a funciton of dimer mean copy number. different behaviors?

from singlegenesystem import *
from SGS import *
import itertools as it
import multiprocessing
from functools import partial
import time


# from matplotlib import gridspec
# import matplotlib.pylab as plt
# if once=True, the code runs trought one random set of parameter and plot stochastic and determinitic time series for the dimer.
# (TO DO) if once=False, the code can be use for large scans, with results saved in txt file (param sets, and time series)
# once = False

def randomscan(r, oscillator, namefiletosavedata, t_f, dt, maxmolecules, shift, system=System.RANDOM):
    # oscillator is the type of model (MDS, 2DS or 3DS),
    # number_of_parameter_sets: size of the scan
    # t_f, dt to define max time of simulation, and saving time.
    savetimeseries = False
    r = r + shift
    # seed for random numbers
    np.random.seed(int(time.time()) + r)

    # parameters that always stay the same, ie no delay.
    pars = {'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'vol': 1}
    # initiate the model
    os = oscillator()
    # rename list of all parameters
    if oscillator == MDS:
        parlist = np.array(['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss',
                   'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm', 'kbd'])
    elif oscillator == DDS:
        parlist = np.array(
            ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12', 'kb1', 'ku1',
             'kb2', 'ku2', 'bcoop', 'ucoop'])
    elif oscillator == DDDS:
        parlist = np.array(
        ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12', 'f13',
         'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop13', 'bcoop23', 'bcoop123',
         'ucoop12', 'ucoop13', 'ucoop23', 'ucoop123'])
    # number of parameters
    numpar = len(parlist)
    # draw random numbers
    # randomparams(oscillator, constraints, Unregulated,SsLrpB=False):
    pars = random_params(oscillator, constraintf=physiological_range, system=system)
    # assign parameters
    os.set_parameters(pars)
    # compute steady-state - the steadystate is stored in variable called "ss"
    os.steady_state()

    # name of files to save parameters leading to multistability or too high steady state (for which we do NOT run stochastic simulations)
    fsteadystatestoohigh = namefiletosavedata + "_too_high.txt"
    fmulti = namefiletosavedata + "_multi.txt"

    # we will only keep track of solution without multistability
    if len(os.ss['m']) == 1:
        # we simulate trajectories only for a reasonable number of molecules
        if system == System.SSLRPB:
            extracondition = os.SsLrpB_compatible()
        else:
            extracondition = True

        if (os.ss['m'][0] < maxmolecules) and (os.ss['d'][0] < maxmolecules) and extracondition:
            print('I will run a simulation')
            print('non monotonicity type', os.monotonicity_type())
            type_monotonicty = os.monotonicity_type()
            # actually doing the stochastic simulation:
            stochtemp = os.stochastic_time_series(t_f, dt, t_i=0, hf={})
            # create empty arrays to save means and variances for each dynamical variable
            meansvars = np.zeros(len(os.allvars))
            varsvars = np.zeros(len(os.allvars))
            fanovars = np.zeros(len(os.allvars))
            # compute mean and variance of each dynamical variable:
            for idx, par in enumerate(os.allvars):
                meansvars[idx] = np.mean(stochtemp.sol[par])
                varsvars[idx] = np.var(stochtemp.sol[par])
            # fanovars[idx]=np.mean(stochtemp.sol[par])/np.var(stochtemp.sol[par])

            #####################################################################################################
            # SAVE DATA
            #####################################################################################################

            # save steady states:
            fsteadystates = namefiletosavedata + "_ss_" + str(type_monotonicty) + "_" + str(r) + ".txt"
            with open(fsteadystates, "a") as file:
                s = ''
                for myvar in os.allvars:
                    s += '%f\t' % (os.ss[myvar][0])
                s += '\n'
                file.write(s)

            # save parameter set considered:
            fpar = namefiletosavedata + "_parms_" + str(type_monotonicty) + "_" + str(r) + ".txt"
            with open(fpar, "a") as file:
                file.write(os.parameter_line())
            file.close()

            # save time series
            if savetimeseries == True:
                ftimeseries = namefiletosavedata + "_ts_" + str(type_monotonicty) + "_" + str(r) + ".txt"
                with open(ftimeseries, "a") as file:
                    s = ''
                    for myvar in os.allvars:
                        # print(myvar)
                        for kk in stochtemp.sol[myvar]:
                            # print(k)
                            s += '%d\t' % kk
                        s += '\n'
                    file.write(s)
                file.close()

            # save mean values of all dynamical variables
            fmeans = namefiletosavedata + "_means_" + str(type_monotonicty) + "_" + str(r) + ".txt"
            with open(fmeans, "a") as file:
                s = ''
                for kk in range(len(meansvars)):
                    temp = meansvars[kk]
                    s += '%f\t' % temp
                s += '\n'
                file.write(s)
            file.close()

            # save variance values of all dynamical variables
            fvars = namefiletosavedata + "_vars_" + str(type_monotonicty) + "_" + str(r) + ".txt"
            with open(fvars, "a") as file:
                s = ''
                for kk in range(len(meansvars)):
                    s += '%f\t' % (varsvars[kk])
                s += '\n'
                file.write(s)
            file.close()

        else:
            # steadystate too high

            with open(fsteadystatestoohigh, "a+") as file:
                s = ''
                for myvar in os.ALL_VARS:
                    s += '%f\t' % (os.ss[myvar])
                s += '\n'
                file.write(s)
            file.close()
            fpar = namefiletosavedata + "_params_too_high.txt"
            with open(fpar, "a+") as file:
                file.write(os.parameter_line())
            file.close()
    else:
        print('mutistability')
        with open(fmulti, "a") as file:
            s = 'hello \n'
            file.write(s)
        file.close()
        # save parameter set considered:
        fpar = namefiletosavedata + "_multi_parms.txt"
        with open(fpar, "a+") as file:
            file.write(os.parameter_line())
        file.close()

    # rename the steadystate
    # hf = os.ss
    # choose initial value for stoch sim (actually useless now but we need to specify something, the function stochacstictimelapse creates the inital condition)
    # y0 = np.array([hf['DNA0'][0],hf['DNA1'][0],hf['DNA2'][0],hf['DNA3'][0],hf['DNA12'][0],hf['DNA23'][0],hf['DNA13'][0],hf['mRNA'][0],hf['m'][0],hf['d'][0]])
    # t=np.arange(0, t_f, dt)
    # dde = odeint(os.eqns, y0, t,args=(1.0,1.0))
    # rseed,oscillator, namefiletosavedata,t_f,dt,maxmolecules


def mainf(oscillator, numb_scanned_sets, namefiletosavedata, t_f, dt, numberofcores, maxmolecules, shift,
          system=System.UNREGULATED):
    np.random.seed(int(time.time()))
    p = multiprocessing.Pool(numberofcores)
    randomscan_assigned = partial(randomscan, oscillator=oscillator, namefiletosavedata=namefiletosavedata, t_f=t_f,
                                  dt=dt, maxmolecules=maxmolecules, shift=shift, system=system)
    r = range(numb_scanned_sets)
    # debug mode
    for kkk in range(200):
        randomscan_assigned(0)
    #p.map(randomscan_assigned, r)
