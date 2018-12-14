from MDS import *
from DDS import *
from DDDS import *

def randomparams(oscillator, constraints, SsLrpB=False,Unregulated=False):
    params = {'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'vol': 1}
    for par in ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0']:
        params[par] = 10 ** (np.random.uniform(np.log10(constraints[par][0]), np.log10(constraints[par][1])))
    if (oscillator == DDS):
        for par in ['kb1', 'kb2']:
            params[par] = 10**(np.random.uniform(np.log10(constraints['kb'][0]), np.log10(constraints['kb'][1])))
        for par in ['ku1', 'ku2']:
            params[par] = 10**(np.random.uniform(np.log10(constraints['ku'][0]), np.log10(constraints['ku'][1])))
        for par in ['ucoop', 'bcoop']:
            params[par] = 10 ** (np.random.uniform(np.log10(constraints['coop'][0]), np.log10(constraints['coop'][1])))
        for par in ['f1', 'f2', 'f12']:
            params[par] = 10**(np.random.uniform(np.log10(constraints['f'][0]), np.log10(constraints['f'][1])))
    elif (oscillator == MDS):
        if Unregulated==True:
            for par in ['kbm', 'kbd', 'kum', 'kud', 'fm', 'fd']:
                params[par] = 0
        else:
            for par in ['kbm', 'kbd']:
                params[par] = 10 ** (np.random.uniform(np.log10(constraints['kb'][0]), np.log10(constraints['kb'][1])))
            for par in ['kum', 'kud']:
                params[par] = 10 ** (np.random.uniform(np.log10(constraints['ku'][0]), np.log10(constraints['ku'][1])))
            for par in ['fm', 'fd']:
                params[par] = 10 ** (np.random.uniform(np.log10(constraints['f'][0]), np.log10(constraints['f'][1])))
    elif (oscillator == DDDS):
        if (SsLrpB):
            ct = 1e-3 / 2.4

            params['Kd1'] = 73.5 * ct
            params['Kd2'] = 0.7 * ct
            params['Kd3'] = 49.1 * ct
            params['kb1'] = 10 ** np.random.uniform(np.log10(max(1e-3, 1e-2 * params['Kd1'])),
                                                    np.log10(min(1e2, 1e3 * params['Kd1'])))
            params['kb2'] = 10 ** np.random.uniform(np.log10(max(1e-3, 1e-2 * params['Kd2'])),
                                                    np.log10(min(1e2, 1e3 * params['Kd2'])))
            params['kb3'] = 10 ** np.random.uniform(np.log10(max(1e-3, 1e-2 * params['Kd3'])),
                                                    np.log10(min(1e2, 1e3 * params['Kd3'])))
            params['bcoop12'] = 10 ** np.random.uniform(np.log10(0.05 * 4.1), np.log10(20))
            params['ucoop12'] = params['bcoop12'] / 4.1
            params['bcoop23'] = 10 ** np.random.uniform(np.log10(0.05 * 2.1), np.log10(20))
            params['ucoop23'] = params['bcoop23'] / 2.1
            params['bcoop13'] = 10 ** np.random.uniform(np.log10(0.05 * 7.9), np.log10(20))
            params['ucoop13'] = params['bcoop13'] / 7.9
            params['bcoop123'] = 10 ** np.random.uniform(np.log10(0.05 * 3.1), np.log10(20))
            params['ucoop123'] = params['bcoop123'] / 3.1

        else:
            for par in ['kb1', 'kb2', 'kb3']:
                params[par] = 10 ** (np.random.uniform(np.log10(constraints['kb'][0]), np.log10(constraints['kb'][1])))
            for par in ['ku1', 'ku2', 'ku3']:
                params[par] = 10 ** (np.random.uniform(np.log10(constraints['ku'][0]), np.log10(constraints['ku'][1])))
            for par in ['bcoop12', 'bcoop23', 'bcoop13', 'bcoop123','ucoop12', 'ucoop23', 'ucoop13', 'ucoop123']:
                params[par] = 10 ** np.random.uniform(np.log10(constraints['coop'][0]), np.log10(constraints['coop'][1]))

        for par in ['f1', 'f2', 'f3', 'f12', 'f13', 'f23', 'f123']:
            params[par] = 10 ** (np.random.uniform(np.log10(constraints['f'][0]), np.log10(constraints['f'][1])))

    return params

def goodrangeparams(oscillator, constraints={}):
    if (oscillator == DDS):
        if (len(constraints)>0):
            params = {'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'vol': 1}
            for par in ['beta', 'gammam', 'gammamRNA', 'Kdim', 'phi0', 'bcoop', 'ucoop']:
                params[par] = 10**(np.random.uniform((np.log10(constraints[par][0])+np.log10(constraints[par][1]))/2, np.log10(constraints[par][1])))
            for par in ['gammad', 'alphaass']:
                params[par] = 10**(np.random.uniform(np.log10(constraints[par][0]), (np.log10(constraints[par][0])+np.log10(constraints[par][1]))/2))
            for par in ['kb1', 'kb2']:
                params[par] = 10**(np.random.uniform(np.log10(constraints['kb'][0]), np.log10(constraints['kb'][1])))
            for par in ['ku1', 'ku2']:
                params[par] = 10**(np.random.uniform(np.log10(constraints['ku'][0]), np.log10(constraints['ku'][1])))
            for par in ['f1', 'f2']:
                params[par] = 10**(np.random.uniform((np.log10(constraints['fa'][0])+np.log10(constraints['fa'][1]))/2, np.log10(constraints['fa'][1]))) #todo fa
            for par in ['f12']:
                params[par] = 10**(np.random.uniform(np.log10(constraints['fr'][0]), (np.log10(constraints['fr'][0])+np.log10(constraints['fr'][1]))/2)) #np.log10(max(params['f1'],params['f2']))))#np.log10(constraints['f'][1]))) #todo fr
        else:
            print("no constraints on random pars")
        return params;

def readparams(i, data, oscillator):
    taum = 0.  # time needed to make a protein that can form a dimer (0.5)
    taumRNA = 0.  # time needed to make mRNA that can be translated (transcription-translation coupling) (10)
    if(oscillator == DDDS):
        # fixed params
        params = {'taum': taum, 'taumRNA': taumRNA, 'DNAtot': 1, 'vol': 1};
        if len(data[i]) == 25:
            parlist = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f123',
                       'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop13', 'bcoop23',
                       'bcoop123', 'ucoop12', 'ucoop13', 'ucoop23', 'ucoop123', 'taum', 'taumRNA'];
        elif len(data[i]) == 30:
            parlist = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0',
                       'f1', 'f2', 'f3', 'f12', 'f13', 'f23', 'f123', 'kb1', 'kb2',
                   'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop13', 'bcoop23', 'bcoop123', 'ucoop12', 'ucoop13',
                   'ucoop23', 'ucoop123', 'taum', 'taumRNA'];
        elif len(data[i]) == 28:
            parlist = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3',
                       'f12', 'f13', 'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop13',
                       'bcoop23', 'bcoop123', 'ucoop12', 'ucoop13', 'ucoop23', 'ucoop123'];
    elif(oscillator == DDS):
        params = {'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'vol': 1};
        if(len(data[i])==16):
            parlist = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12',
                       'kb1', 'ku1', 'kb2', 'ku2', 'bcoop', 'ucoop'];
        elif(len(data[i])==17): #17th element is selection
            parlist = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12',
                       'kb1', 'ku1', 'kb2', 'ku2', 'bcoop', 'ucoop'];
        else:
            parlist = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1',  'f12',
                       'kb1', 'ku1', 'kb2', 'ku2', 'bcoop', 'ucoop'];
    elif(oscillator == MDS):
        params = {'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'vol': 1};
        parlist = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm',
               'kbd'];
    elif(oscillator == DS):
        params = {'taum': 0, 'taumRNA': 0, 'DNAtot': 1, 'vol': 1};
        parlist = ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'kb1', 'ku1'];

    # read parameters from file
    for j in range(len(parlist)):
        params[parlist[j]] = data[i][j];
    if(oscillator == DDDS and len(data[i]) == 25):
        params['f2'] = params['f1']
        params['f3'] = params['f1']
        params['f12'] = params['f1']
        params['f23'] = params['f1']
        params['f13'] = params['f1']
    elif(oscillator==DDS and len(data[i]) == 15):
        params['f2'] = params['f1']
    #if(oscillator == DDO):
    #    params['f2'] = params['f12']
    return params