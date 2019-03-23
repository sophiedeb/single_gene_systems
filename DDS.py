#!/usr/bin/env python

"""DDS.py: Class two dimer system (2DS)."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from singlegenesystem import *

class DDS(SG):
    PARAMETERS = np.array(
        ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12', 'kb1', 'ku1',
         'kb2', 'ku2', 'bcoop', 'ucoop', 'DNAtot', 'taum', 'taumRNA'])
    DNA_STATES = np.array(['DNA0', 'DNA1', 'DNA2', 'DNA12'])
    ALL_VARS = np.array(['DNA0', 'DNA1', 'DNA2', 'mRNA', 'm', 'd'])

    def __init__(self, params={}):


        self.ss = {}

        self.set_parameters(params)

    #--------------------------------------------------------------------------------------------------------
    # PARAMETER FUNCTIONS
    #--------------------------------------------------------------------------------------------------------

    def set_parameters(self, params):
        """ Set parameters of the 2DS.

        Only parameters of parameterlist are accepted, other parameters are initialized by None.
        Automatic conversion from binding constants to unbinding rates and from cooperativity to unbinding cooperativity.
        The steady state is cleared, because different parameters will lead to a different steady state.
        """

        for param in self.PARAMETERS:
            if param in params:
                setattr(self, param, params[param])
            else:
                setattr(self, param, None)

        for i in [1, 2]:
            if 'Kd%d' % i in params:
                setattr(self, 'Kd%d' % i, params['Kd%d' % i])
                setattr(self, 'ku%d' % i, params['kb%d' % i] / params['Kd%d' % i])
            elif 'ku%d' % i in params:
                setattr(self, 'ku%d' % i, params['ku%d' % i])
                setattr(self, 'Kd%d' % i, params['kb%d' % i] / params['ku%d' % i])

        if 'omega' in params:
            setattr(self, 'ucoop', params['omega'] / params['bcoop'])

        self.ss.clear()

        #if ('alphadiss' in params):
        #    self.alphadiss = params['alphadiss'];
        #    self.Kdim = self.alphadiss / self.alphaass;
        #elif ('Kdim' in params):
        #    self.Kdim = params['Kdim'];
        #    self.alphadiss = self.alphaass * self.Kdim;

    def __copy__(self):
        """ Makes a copy."""

        newone = type(self)()
        newone.set_parameters(self.get_parameters())

        return newone

    # --------------------------------------------------------------------------------------------------------
    # EQUATIONS AND TIME SERIES
    # --------------------------------------------------------------------------------------------------------

    def eqns(self, var, t=0.0, fastdimer=False):
        """ Ordinary differential equations of system.
        """

        if not fastdimer:
            # Set variables
            DNA0, DNA1, DNA2, mRNA, m, d = var

            # Calculate all differential equations.
            dDNA0 = (self.ku1 * DNA1 + self.ku2 * DNA2) - (self.kb1 + self.kb2) * d * DNA0
            dDNA1 = self.kb1 * d * DNA0 + self.ucoop * self.ku2 * (self.DNAtot - DNA0 - DNA1 - DNA2) - (
                                                                                                           self.ku1 + self.bcoop * self.kb2 * d) * DNA1
            dDNA2 = self.kb2 * d * DNA0 + self.ucoop * self.ku1 * (self.DNAtot - DNA0 - DNA1 - DNA2) - (
                                                                                                           self.ku2 + self.bcoop * self.kb1 * d) * DNA2
            dmRNA = self.phi0 * ((1 - self.f12)*DNA0 + (self.f1 - self.f12) * DNA1 + (self.f2 - self.f12) * DNA2 + self.f12 * self.DNAtot) - self.gammamRNA * mRNA
            dm = self.beta * mRNA - 2 * self.alphaass * pow(m, 2) + 2 * self.alphadiss * d - self.gammam * m
            dd = -d * ((self.kb1 + self.kb2) * DNA0 + self.bcoop * (self.kb2 * DNA1 + self.kb1 * DNA2)) + (
                self.ku1 * DNA1 + self.ku2 * DNA2 + (self.ku1 + self.ku2) * self.ucoop * (
                    self.DNAtot - DNA0 - DNA1 - DNA2)) + self.alphaass * pow(m, 2) - self.alphadiss * d - self.gammad * d
            return np.array([dDNA0, dDNA1, dDNA2, dmRNA, dm, dd])
        else:
            # Set variables.
            mRNA, m = var

            qss = self.quasi_steady_state(m)
            DNA0 = qss['DNA0']
            DNA1 = qss['DNA1']
            DNA2 = qss['DNA2']
            d = qss['d']

            dmRNA = self.phi0 * ((1 - self.f12)*DNA0 + (self.f1 - self.f12) * DNA1 + (self.f2 - self.f12) * DNA2 + self.f12 * self.DNAtot) - self.gammamRNA * mRNA
            dm = self.beta * mRNA - 2 * self.alphaass * pow(m, 2) + 2 * self.alphadiss * d - self.gammam * m
            return [dmRNA, dm]

    @staticmethod
    def reaction_matrix():
        """ Returns the reaction matrix.

        Columns are chemical processes. Rows are the change of every variable at given process.
        """

        # Columns are chemical processes. For order see function propensities.
        # Rows are the change of the variables after specific chemical process.
        # Rows in order : DNA0,
        reactions = np.matrix(
            [[-1, 1,  -1, 1,  0,  0,  0,  0,  0, 0, 0,  0,  0,  0,  0],
             [1,  -1, 0,  0,  -1, 1,  0,  0,  0, 0, 0,  0,  0,  0,  0],
             [0,  0,  1,  -1, 0,  0,  -1, 1,  0, 0, 0,  0,  0,  0,  0],
             [0,  0,  0,  0,  0,  0,  0,  0,  1, 0, 0,  0,  -1, 0,  0],
             [0,  0,  0,  0,  0,  0,  0,  0,  0, 1, -2, 2,  0,  -1, 0],
             [-1, 1,  -1, 1,  -1, 1,  -1, 1,  0, 0, 1,  -1, 0,  0,  -1]]
        )
        return reactions

    def propensities(self, var):
        """ Propensities for variables = arg."""

        # Get variables.
        DNA0, DNA1, DNA2, mRNA, m, d = var

        # Set propensities for all processes.
        prop = np.array(
                  [self.kb1*DNA0*d, # DNA0 + d --> DNA1 (binding)
                  self.ku1*DNA1, # DNA1 --> DNA0 + d (unbinding)
                  self.kb2*DNA0*d, # DNA0 + d --> DNA2 (binding)
                  self.ku2*DNA2, # DNA2 --> DNA0 + d (unbinding)
                  self.kb2*self.bcoop*DNA1*d, # DNA1 + d --> DNA12 (binding)
                  self.ku2*self.ucoop*(self.DNAtot-DNA0-DNA1-DNA2), # DNA12 --> DNA1 + d (unbinding)
                  self.kb1*self.bcoop*DNA2*d, # DNA2 + d --> DNA12 (binding)
                  self.ku1*self.ucoop*(self.DNAtot-DNA0-DNA1-DNA2), # DNA12 --> DNA2 + d (unbinding)
                  self.phi0*((1 - self.f12)*DNA0 + (self.f1 - self.f12)*DNA1
                             + (self.f2 - self.f12)*DNA2 + self.f12*self.DNAtot), # --> mRNA (transcription)
                  self.beta*mRNA, # mRNA --> mRNA + m (translation)
                  self.alphaass*m*(m-1), # 2*m --> d (dimer association)
                  self.alphadiss*d, # d --> 2*m (dimer dissociation)
                  self.gammamRNA*mRNA, # mRNA --> 0 (mRNA degradation)
                  self.gammam*m, # m --> 0 (monomer degradation)
                  self.gammad*d]) # d --> 0 (dimer degradation)

        # Print message if any propensity is negative.
        if np.any(prop<0):
            print(prop)

        return prop

    # --------------------------------------------------------------------------------------------------------
    # STEADY STATE, JACOBIAN & EIGENVALUES
    # --------------------------------------------------------------------------------------------------------

    def set_steady_state(self):
        """ Set the steady state of the system."""

        # Define all auxiliary functions.
        fd = lambda x: self.alphaass * x ** 2 / (self.alphadiss + self.gammad)
        f0 = lambda x: self.DNAtot / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)
        f1 = lambda x: self.DNAtot * fd(x) * self.Kd1 / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)
        f2 = lambda x: self.DNAtot * fd(x) * self.Kd2 / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)

        func = lambda x: (self.phi0 * ((1 - self.f12) * f0(x) + (self.f1 - self.f12) * f1(x)
                                       + (self.f2 - self.f12) * f2(x) + self.f12 * self.DNAtot) - self.gammamRNA / self.beta * (
                              self.gammam * x) - self.gammamRNA / self.beta * 2 * (
                              self.alphaass * x ** 2 - self.alphadiss * fd(x)))

        # Find steady state(s) of the monomer.
        m = self._find_x_steady_state(func)

        # Calculate the steady state values for every variable from the monomer steady state.
        d = fd(m)
        DNA0 = f0(m)
        DNA1 = f1(m)
        DNA2 = f2(m)
        mRNA = self.phi0 * ((1 - self.f12)*DNA0 + (self.f1 - self.f12) * DNA1 + (self.f2 - self.f12) * DNA2 + self.f12 * self.DNAtot) / self.gammamRNA

        # Set all variables of the steady state.
        self.ss = {'DNA0': DNA0, 'DNA1': DNA1, 'DNA2': DNA2, 'DNA12': self.DNAtot - DNA0 - DNA1 - DNA2, 'mRNA': mRNA,
                   'm': m, 'd': d}

    def quasi_steady_state(self, x, polymer=Polymer.MONOMER):
        """ Returns a dictionary with the quasi steady state for x, the monomer or polymer count.

                steady state --> assuming fast (un)binding and dimerization:
                solution to dm/dt = 0 and dmRNA/dt = 0 under the assumption that dDNAi/dt and dd/dt = 0.
                """

        # Convert monomer and dimer rates.
        if polymer == Polymer.MONOMER:
            m = x
            d = self.alphaass * m ** 2 / (self.alphadiss + self.gammad)
        elif polymer == Polymer.DIMER:
            d = x
            m = np.sqrt((self.alphadiss + self.gammad) * d / self.alphaass)

        # Calculate DNA states given the monomer concentration.
        fd = lambda x: self.alphaass * x ** 2 / (self.alphadiss + self.gammad)
        f0 = lambda x: self.DNAtot / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)
        f1 = lambda x: self.DNAtot * fd(x) * self.Kd1 / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)
        f2 = lambda x: self.DNAtot * fd(x) * self.Kd2 / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)

        DNA0 = f0(m)
        DNA1 = f1(m)
        DNA2 = f2(m)
        qss = {'DNA0': DNA0, 'DNA1': DNA1, 'DNA2': DNA2, 'DNA12': self.DNAtot - DNA0 - DNA1 - DNA2, 'd': d, 'm': m}

        return qss

    def jacobian(self, var):
        # m, d, DNA0, DNA1, DNA2 = var;
        """ Return the jacobian of the system for variables = arg. """

        # Get the variables.
        DNA0, DNA1, DNA2, mRNA, m, d = var

        # Define the Jacobian.
        J = [[-(self.kb1 + self.kb2) * d, self.ku1, self.ku2, 0, 0, -(self.kb1 + self.kb2) * DNA0],
             [self.kb1 * d - self.ucoop * self.ku2, -self.ucoop * self.ku2 - (self.ku1 + self.bcoop * self.kb2 * d),
              -self.ucoop * self.ku2, 0, 0, self.kb1 * DNA0 - self.bcoop * self.kb2 * DNA1],
             [self.kb2 * d - self.ucoop * self.ku1, -self.ucoop * self.ku1,
              -self.ucoop * self.ku1 - (self.ku2 + self.bcoop * self.kb1 * d), 0, 0,
              self.kb2 * DNA0 - self.bcoop * self.kb1 * DNA2],
             [self.phi0 * (1 - self.f12), self.phi0 * (self.f1 - self.f12), self.phi0 * (self.f2 - self.f12), -self.gammamRNA, 0, 0],
             [0, 0, 0, self.beta, -4 * self.alphaass * m - self.gammam, 2 * self.alphadiss],
             [-(self.kb1 + self.kb2) * d - (self.ku1 + self.ku2) * self.ucoop,
              -self.kb2 * d * self.bcoop + (self.ku1 - (self.ku1 + self.ku2) * self.ucoop),
              -self.kb1 * d * self.bcoop + (self.ku2 - (self.ku1 + self.ku2) * self.ucoop), 0, 2 * self.alphaass * m,
              -self.alphadiss - self.gammad - (
                  (self.kb1 + self.kb2) * DNA0 + self.kb2 * self.bcoop * DNA1 + self.kb1 * self.bcoop * DNA2)]]

        return J

    # --------------------------------------------------------------------------------------------------------
    # AUXILIARY FUNCTIONS
    # --------------------------------------------------------------------------------------------------------

    """
    def transfunc(self):
        a = (self.Kd1*self.Kd2*self.bcoop/self.ucoop*(self.Kd1*(self.f12 - self.f1)+ self.Kd2*(self.f12 - self.f2)))
        b = (2*self.Kd1*self.Kd2*self.bcoop/self.ucoop*(self.f12 - 1))
        c = self.Kd1*(self.f1 - 1) + self.Kd2*(self.f2 - 1)
        D = b**2 - 4*a*c
        x1 = (-b-np.sqrt(D))/(2*a) if D>0 else np.nan
        x2 = (-b+np.sqrt(D))/(2*a) if D>0 else np.nan
        return D, x1, x2
    """

    def dimensionless_parameters(self, ts, qs):
        """ Return parameters without dimensions using time scale ts and quantity scale qs."""

        par = {  # nondimensionalized
            'beta': self.beta * ts, 'gammam': self.gammam * ts, 'gammamRNA': self.gammamRNA * ts,
            'gammad': self.gammad * ts, 'alphaass': self.alphaass * ts * qs, 'alphadiss': self.alphadiss * ts,
            'taum': self.taum / ts, 'taumRNA': self.taumRNA / ts, 'DNAtot': self.DNAtot / qs, 'phi0': self.phi0 * ts,
            'fa': self.fa, 'fr': self.fr, 'Km': self.Km * qs, 'Kd': self.Kd * qs, 'kbm': self.kbm * qs * ts,
            'kum': self.kum * ts, 'kb': self.kb * qs * ts, 'ku': self.ku * ts, 'volume': self.vol, 'ts': ts, 'qs': qs}

        return par

    """
    def sscurve(self, ax, e=5000, c='k'):
        if (not len(self.ss) > 0):
            self.set_steady_state()

        fd = lambda x: self.alphaass * x ** 2 / (self.alphadiss + self.gammad)
        f0 = lambda x: self.DNAtot / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)
        f1 = lambda x: self.DNAtot * fd(x) * self.Kd1 / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)
        f2 = lambda x: self.DNAtot * fd(x) * self.Kd2 / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)

        func = lambda x: (self.phi0 * (
        (1 - self.f12) * f0(x) + (self.f1 - self.f12) * f1(x) + (self.f2 - self.f12) * f2(
            x) + self.f12 * self.DNAtot) - self.gammamRNA / self.beta * (
                              self.gammam * x) - self.gammamRNA / self.beta * 2 * (
                              self.alphaass * x ** 2 - self.alphadiss * fd(x)))

        x = np.linspace(0, e, 2000)
        fx = func(x)
        if(len(self.ss['d'])>1):
            posx = fx[np.logical_and(x > self.ss['d'][1], x < self.ss['d'][2])]
            score = sum(fx)/len(fx)
            ax.plot(fd(x), fx, color=c, label=('%.2f'%score))
        #ax.set_xscale('log')
        #ax.set_ylim([-10, 10])


        ax.axhline(0, color='k', linestyle=':')
    """