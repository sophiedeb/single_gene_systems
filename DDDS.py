#!/usr/bin/env python

"""DDDS.py: Class three dimer system (3DS)."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from singlegenesystem import *

class DDDS(SG):
    PARAMETERS = np.array(
        ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12', 'f13',
         'f23', 'f123', 'kb1', 'kb2', 'kb3', 'ku1', 'ku2', 'ku3', 'bcoop12', 'bcoop13', 'bcoop23', 'bcoop123',
         'ucoop12', 'ucoop13', 'ucoop23', 'ucoop123', 'DNAtot', 'taum', 'taumRNA'])

    DNA_STATES = ['DNA0', 'DNA1', 'DNA2', 'DNA3', 'DNA12', 'DNA23', 'DNA13', 'DNA123']
    ALL_VARS = np.array(['DNA0', 'DNA1', 'DNA2', 'DNA3', 'DNA12', 'DNA23', 'DNA13', 'mRNA', 'm', 'd'])

    def __init__(self, params={}):
        self.ss = {}  # steady state

        self.set_parameters(params)

    #--------------------------------------------------------------------------------------------------------
    # PARAMETER FUNCTIONS
    #--------------------------------------------------------------------------------------------------------

    def set_parameters(self, params):
        """ Set parameters of the 3DS.

        Only parameters of parameterlist are accepted, other parameters are initialized by None.
        Automatic conversion from binding constants to unbinding rates and from cooperativity to unbinding cooperativity.
        The steady state is cleared, because different parameters will lead to a different steady state.
        """

        for param in self.PARAMETERS:
            if param in params:
                setattr(self, param, params[param])
            else:
                setattr(self, param, None)

        for i in range(1, 4):
            if 'Kd%d'%i in params:
                setattr(self, 'ku%d'%i, params['kb%d'%i] / params['Kd%d'%i])

        for i in [12, 23, 13, 123]:
            if 'omega%d'%i in params:
                setattr(self, 'ucoop%d'%i, params['omega%d'%i] / params['bcoop%d'%i])

        self.ss.clear()

    #def __copy__(self):
    #    """ Makes a copy."""

    #        newone = type(self)()
    #    newone.set_parameters(self.get_parameters())

    #    return newone

    def small_mutant(self, parent, constraintf, withdelay=False, system=System.RANDOM):
        """ Set parameters to the one of parent with one value slightly changed within the constraints by the constraintfunction.
        """

        # Small mutant as any single gene system.
        super().small_mutant(parent, constraintf, withdelay)

        # For SsLrpB parameters repair the known parameters.
        if system == System.SSLRPB:
            self.ku1 = self.kb1 / 1.764
            self.ku2 = self.kb2 / 0.0168
            self.ku3 = self.kb3 / 1.1784
            self.ucoop12 = self.bcoop12 / 4.1
            self.ucoop23 = self.bcoop23 / 7.9
            self.ucoop13 = self.bcoop13 / 2.1
            self.ucoop123 = self.bcoop123 / 3.1

        return

    def large_mutant(self, parent, constraintf, withdelay=False, system=System.RANDOM):
        """ Set parameters to the one of parent with one value randomly changed within the constraints by the constraintfunction.
        """

        # Large mutant as any single gene system.
        super().large_mutant(parent, constraintf, withdelay)

        # For SsLrpB parameters repair the known parameters.
        if system == System.SSLRPB:
            self.ku1 = self.kb1 / 1.764
            self.ku2 = self.kb2 / 0.0168
            self.ku3 = self.kb3 / 1.1784
            self.ucoop12 = self.bcoop12 / 4.1
            self.ucoop23 = self.bcoop23 / 7.9
            self.ucoop13 = self.bcoop13 / 2.1
            self.ucoop123 = self.bcoop123 / 3.1

        return

    # --------------------------------------------------------------------------------------------------------
    # EQUATIONS AND TIME SERIES
    # --------------------------------------------------------------------------------------------------------

    def eqns(self, var, t):
        """ Ordinary differential equations of system.
        """

        # Set variables
        DNA0, DNA1, DNA2, DNA3, DNA12, DNA23, DNA13, mRNA, m, d = var

        # Calculate all differential equations.
        dDNA0 = (self.ku1 * DNA1 + self.ku2 * DNA2 + self.ku3 * DNA3) - (self.kb1 + self.kb2 + self.kb3) * d * DNA0
        dDNA1 = self.kb1 * d * DNA0 + self.ucoop12 * self.ku2 * DNA12 + self.ucoop13 * self.ku3 * DNA13 - (self.ku1 + (
            self.bcoop12 * self.kb2 + self.bcoop13 * self.kb3) * d) * DNA1
        dDNA2 = self.kb2 * d * DNA0 + self.ucoop12 * self.ku1 * DNA12 + self.ucoop23 * self.ku3 * DNA23 - (self.ku2 + (
            self.bcoop12 * self.kb1 + self.bcoop23 * self.kb3) * d) * DNA2
        dDNA3 = self.kb3 * d * DNA0 + self.ucoop13 * self.ku1 * DNA13 + self.ucoop23 * self.ku2 * DNA23 - (self.ku3 + (
            self.bcoop13 * self.kb1 + self.bcoop23 * self.kb2) * d) * DNA3
        dDNA12 = self.bcoop12 * d * (
            self.kb1 * DNA2 + self.kb2 * DNA1) + self.ucoop123 * self.ucoop13 * self.ucoop23 * self.ku3 * (
            self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (
                                                                                   self.bcoop123 * self.bcoop13 * self.bcoop23 * self.kb3 * d + self.ucoop12 * (
                                                                                       self.ku1 + self.ku2)) * DNA12
        dDNA23 = self.bcoop23 * d * (
            self.kb3 * DNA2 + self.kb2 * DNA3) + self.ucoop123 * self.ucoop12 * self.ucoop13 * self.ku1 * (
            self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (
                                                                                   self.bcoop123 * self.bcoop12 * self.bcoop13 * self.kb1 * d + self.ucoop23 * (
                                                                                       self.ku2 + self.ku3)) * DNA23
        dDNA13 = self.bcoop13 * d * (
            self.kb1 * DNA3 + self.kb3 * DNA1) + self.ucoop123 * self.ucoop12 * self.ucoop23 * self.ku2 * (
            self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (
                                                                                   self.bcoop123 * self.bcoop12 * self.bcoop23 * self.kb2 * d + self.ucoop13 * (
                                                                                       self.ku1 + self.ku3)) * DNA13
        dmRNA = self.phi0 * (
        DNA0 + self.f1 * DNA1 + self.f2 * DNA2 + self.f3 * DNA3 + self.f12 * DNA12 + self.f23 * DNA23 + self.f13 * DNA13 + self.f123 * (
            self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23)) - self.gammamRNA * mRNA
        dm = self.beta * mRNA - 2 * self.alphaass * pow(m, 2) + 2 * self.alphadiss * d - self.gammam * m
        dd = (-d * (
            (self.kb1 + self.kb2 + self.kb3) * DNA0 + (self.kb2 * self.bcoop12 + self.kb3 * self.bcoop13) * DNA1 + (
                self.kb1 * self.bcoop12 + self.kb3 * self.bcoop23) * DNA2 + (
                self.kb1 * self.bcoop13 + self.kb2 * self.bcoop23) * DNA3 + self.kb3 * self.bcoop123 * self.bcoop13 * self.bcoop23 * DNA12 + self.kb2 * self.bcoop123 * self.bcoop12 * self.bcoop23 * DNA13 + self.kb1 * self.bcoop123 * self.bcoop12 * self.bcoop13 * DNA23) + (
                  self.ku1 * DNA1 + self.ku2 * DNA2 + self.ku3 * DNA3 + (self.ku1 + self.ku2) * self.ucoop12 * DNA12 + (
                  self.ku1 + self.ku3) * self.ucoop13 * DNA13 + (
                      self.ku2 + self.ku3) * self.ucoop23 * DNA23 + self.ucoop123 * (
                      self.ku1 * self.ucoop12 * self.ucoop13 + self.ku2 * self.ucoop12 * self.ucoop23 + self.ku3 * self.ucoop13 * self.ucoop23) * (
                      self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23)) + self.alphaass * pow(m,
                                                                                                              2) - self.alphadiss * d - self.gammad * d)
        return [dDNA0, dDNA1, dDNA2, dDNA3, dDNA12, dDNA23, dDNA13, dmRNA, dm, dd]

    @staticmethod
    def eqns_delay():
        """ Delay differential equations.

        Use of a dictionary of functions for package pydelay.
        """

        eq = {'DNA0': 'ku1 * DNA1 + ku2 * DNA2 + ku3 * DNA3 - (kb1 + kb2 + kb3) * d * DNA0',
              'DNA1': 'kb1 * d * DNA0 + ucoop12 * ku2 * DNA12 + ucoop13 * ku3 * DNA13 - (ku1 + (bcoop12 * kb2 + bcoop13 * kb3) * d) * DNA1',
              'DNA2': 'kb2 * d * DNA0 + ucoop12 * ku1 * DNA12 + ucoop23 * ku3 * DNA23 - (ku2 + (bcoop12 * kb1 + bcoop23 * kb3) * d) * DNA2',
              'DNA3': 'kb3 * d * DNA0 + ucoop13 * ku1 * DNA13 + ucoop23 * ku2 * DNA23 - (ku3 + (bcoop13 * kb1 + bcoop23 * kb2) * d) * DNA3',
              'DNA12': 'bcoop12 * d * (kb1 * DNA2 + kb2 * DNA1) + ucoop123 * ucoop13 * ucoop23 * ku3 * (DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (bcoop123 * bcoop13 * bcoop23 * kb3 * d + ucoop12 * (ku1 + ku2)) * DNA12',
              'DNA13': 'bcoop13 * d * (kb1 * DNA3 + kb3 * DNA1) + ucoop123 * ucoop12 * ucoop23 * ku2 * (DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (bcoop123 * bcoop12 * bcoop23 * kb2 * d + ucoop13 * (ku1 + ku3)) * DNA13',
              'DNA23': 'bcoop23 * d * (kb3 * DNA2 + kb2 * DNA3) + ucoop123 * ucoop12 * ucoop13 * ku1 * (DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (bcoop123 * bcoop12 * bcoop13 * kb1 * d + ucoop23 * (ku2 + ku3)) * DNA23',
              'mRNA': 'phi0 * (DNA0 + fa * (DNA1 + DNA2 + DNA3 + DNA12 + DNA23 + DNA13) + fr * (DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23)) - gammamRNA * mRNA',
              'm': 'beta * mRNA - 2 * alphaass * pow(m, 2) + 2 * alphadiss * d - gammam * m',
              'd': '(-d * ((kb1 + kb2 + kb3) * DNA0 + (kb2 * bcoop12 + kb3 * bcoop13) * DNA1 + (kb1 * bcoop12 + kb3 * bcoop23) * DNA2 + (kb1 * bcoop13 + kb2 * bcoop23) * DNA3 + kb3 * bcoop123 * bcoop13 * bcoop23 * DNA12 + kb2 * bcoop123 * bcoop12 * bcoop23 * DNA13 + kb1 * bcoop123 * bcoop12 * bcoop13 * DNA23) + (ku1 * DNA1 + ku2 * DNA2 + ku3 * DNA3 + (ku1 + ku2) * ucoop12 * DNA12 + (ku1 + ku3) * ucoop13 * DNA13 + (ku2 + ku3) * ucoop23 * DNA23 + ucoop123 * (ku1 * ucoop12 * ucoop13 + ku2 * ucoop12 * ucoop23 + ku3 * ucoop13 * ucoop23) * (DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23)) + alphaass * pow(m,2) - alphadiss * d - gammad * d)'}
        return eq

    @staticmethod
    def reaction_matrix():
        """ Returns the reaction matrix.

        Columns are chemical processes. Rows are the change of every variable at given process.
        """

        # Columns are chemical processes. For order see function propensities.
        # Rows are the change of the variables after specific chemical process.
        # Rows in order : DNA0,
        reactions = np.matrix(
            [[-1,   1, -1,  1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
             [1,   -1,  0,  0,  0,  0, -1,  1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
             [0,    0,  1, -1,  0,  0,  0,  0,  0,  0, -1,  1,  -1, 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
             [0,    0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0, -1,  1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
             [0,    0,  0,  0,  0,  0,  1, -1,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
             [0,    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  -1, 0,  0,  1, -1,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
             [0,    0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  0,  0],
             [0,    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0, -1,  0,  0],
             [0,    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  2,  0, -1,  0],
             [-1,   1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1,  -1, 1, -1,  1, -1,  1, -1,  1, -1,  1,  -1, 1,  0,  0,  1, -1,  0,  0, -1]]
        )
        return reactions

    def propensities(self, var):
        """ Propensities for variables = arg."""

        # Get variables.
        DNA0, DNA1, DNA2, DNA3, DNA12, DNA23, DNA13, mRNA, m, d = var
        DNA123 = self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13

        # Set propensities for all processes.
        prop = np.array(
                  [self.kb1*DNA0*d, # DNA0 + d --> DNA1 (binding)
                    self.ku1*DNA1, # DNA1 --> DNA0 + d (unbinding)
                    self.kb2*DNA0*d, # DNA0 + d --> DNA2 (binding)
                    self.ku2*DNA2, # DNA2 --> DNA0 + d (unbinding)
                    self.kb3 * DNA0 * d, # DNA0 + d --> DNA3 (binding)
                    self.ku3 * DNA3, # DNA3 --> DNA0 + d (unbinding)
                    self.kb2*self.bcoop12*DNA1*d, # DNA1 + d --> DNA12 (binding)
                    self.ku2*self.ucoop12*DNA12, # DNA12 --> DNA1 + d (unbinding)
                    self.kb3 * self.bcoop13 * DNA1 * d, # DNA1 + d --> DNA13 (binding)
                    self.ku3 * self.ucoop13 * DNA13, # DNA13 --> DNA1 + d (unbinding)
                    self.kb1 * self.bcoop12 * DNA2 * d, # DNA2 + d --> DNA12 (binding)
                    self.ku1 * self.ucoop12 * DNA12, # DNA12 --> DNA2 + d (unbinding)
                    self.kb3 * self.bcoop23 * DNA2 * d, # DNA2 + d --> DNA23 (binding)
                    self.ku3 * self.ucoop23 * DNA23, # DNA23 --> DNA2 + d (unbinding)
                    self.kb1 * self.bcoop13 * DNA3 * d, # DNA3 + d --> DNA13 (binding)
                    self.ku1 * self.ucoop13 * DNA13, # DNA13 --> DNA3 + d (unbinding)
                    self.kb2 * self.bcoop23 * DNA3 * d, # DNA3 + d --> DNA23 (binding)
                    self.ku2 * self.ucoop23 * DNA23, # DNA23 --> DNA3 + d (unbinding)
                    self.kb3 * self.bcoop13 * self.bcoop23 * self.bcoop123 * DNA12*d, # DNA12 + d --> DNA123 (binding)
                    self.ku3 * self.ucoop13 * self.ucoop23 * self.ucoop123 * DNA123, # DNA123 --> DNA12 + d (unbinding)
                    self.kb1 * self.bcoop13 * self.bcoop12 * self.bcoop123 * DNA23 * d, # DNA23 + d --> DNA123 (binding)
                    self.ku1 * self.ucoop13 * self.ucoop12 * self.ucoop123 * DNA123, # DNA123 --> DNA23 + d (unbinding)
                    self.kb2 * self.bcoop12 * self.bcoop23 * self.bcoop123 * DNA13 * d, # DNA13 + d --> DNA123 (binding)
                    self.ku2 * self.ucoop12 * self.ucoop23 * self.ucoop123 * DNA123, # DNA123 --> DNA13 + d (unbinding)
                    self.phi0*(DNA0 + self.f1*DNA1 + self.f2*DNA2 + self.f3*DNA3
                               + self.f12*DNA12 + self.f13*DNA13 + self.f23*DNA23
                               + self.f123*DNA123), # --> mRNA (transcription)
                    self.beta*mRNA, # mRNA --> mRNA + m (translation)
                    self.alphaass*m*(m-1), # 2*m --> d (dimer association)
                    self.alphadiss*d, # d --> 2*m (dimer dissociation)
                    self.gammamRNA*mRNA, # mRNA --> 0 (mRNA degradation)
                    self.gammam*m, # m --> 0 (monomer degradation)
                    self.gammad*d]) # d --> 0 (dimer degradation)

        # Print message if any propensity is negative.
        if np.any(prop<0):
            print("Negative propensity:", prop)

        return prop

    # --------------------------------------------------------------------------------------------------------
    # STEADY STATE, JACOBIAN & EIGENVALUES
    # --------------------------------------------------------------------------------------------------------

    def set_steady_state(self):
        """ Set the steady state of the system."""

        # Define all auxiliary functions.
        fdenom = lambda d: (self.kb1 * self.kb2 * self.kb3 / (
            self.ku1 * self.ku2 * self.ku3) * self.bcoop12 * self.bcoop23 * self.bcoop13 * self.bcoop123 / (
                                self.ucoop12 * self.ucoop23 * self.ucoop13 * self.ucoop123) * d ** 3 + (
                                self.kb1 * self.kb2 / (
                                    self.ku1 * self.ku2) * self.bcoop12 / self.ucoop12 + self.kb1 * self.kb3 / (
                                    self.ku1 * self.ku3) * self.bcoop13 / self.ucoop13 + self.kb2 * self.kb3 / (
                                    self.ku2 * self.ku3) * self.bcoop23 / self.ucoop23) * d ** 2 + (
                                self.kb1 / self.ku1 + self.kb2 / self.ku2 + self.kb3 / self.ku3) * d + 1)
        fd = lambda m: self.alphaass * m ** 2 / (self.alphadiss + self.gammad)
        fDNA0 = lambda d: self.DNAtot / fdenom(d)
        fDNA1 = lambda d: self.DNAtot * self.kb1 / self.ku1 * d / fdenom(d)
        fDNA2 = lambda d: self.DNAtot * self.kb2 / self.ku2 * d / fdenom(d)
        fDNA3 = lambda d: self.DNAtot * self.kb3 / self.ku3 * d / fdenom(d)
        fDNA12 = lambda d: self.DNAtot * self.kb1 * self.kb2 / (
            self.ku1 * self.ku2) * self.bcoop12 / self.ucoop12 * d ** 2 / fdenom(d)
        fDNA13 = lambda d: self.DNAtot * self.kb1 * self.kb3 / (
            self.ku1 * self.ku3) * self.bcoop13 / self.ucoop13 * d ** 2 / fdenom(d)
        fDNA23 = lambda d: self.DNAtot * self.kb2 * self.kb3 / (
            self.ku2 * self.ku3) * self.bcoop23 / self.ucoop23 * d ** 2 / fdenom(d)

        # Function of the monomer concentration for dmRNA/dt = 0.
        func = lambda x: (self.phi0 * ((1 - self.f123)*fDNA0(fd(x)) + (self.f1 - self.f123) * fDNA1(fd(x))
                                       + (self.f2 - self.f123)*fDNA2(fd(x)) + (self.f3 - self.f123)*fDNA3(fd(x))
                        + (self.f12 - self.f123)*fDNA12(fd(x)) + (self.f13 - self.f123)*fDNA13(fd(x))
                        + (self.f23 - self.f123)*fDNA23(fd(x)) + self.f123 * self.DNAtot)
                          - self.gammamRNA / self.beta * (
                              self.gammam * x) - self.gammamRNA / self.beta * 2 * self.alphaass * x ** 2 / (
                              self.alphadiss / self.gammad + 1))

        # Find steady state(s) of the monomer.
        m = self._find_x_steady_state(func)

        # Control plot to see if right steady state value was found.
        if False:
            x = np.logspace(-19, 25, 500)
            plt.figure()
            plt.plot(x, func(x))
            plt.xscale('log')
            plt.ylim([-0.1, 0.1])
            plt.plot(m, 0, '*')
            plt.show()

        # Calculate the steady state values for every variable from the monomer steady state.
        d = fd(m)
        DNA0 = fDNA0(d)
        DNA1 = fDNA1(d)
        DNA2 = fDNA2(d)
        DNA3 = fDNA3(d)
        DNA12 = fDNA12(d)
        DNA13 = fDNA13(d)
        DNA23 = fDNA23(d)
        DNA123 = self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23
        mRNA = self.phi0 * ((1 - self.f123) * DNA0 + (self.f1 - self.f123) * DNA1
                            + (self.f2 - self.f123) * DNA2 + (self.f3 - self.f123) * DNA3
                            + (self.f12 - self.f123) * DNA12 + (self.f13 - self.f123) * DNA13
                            + (self.f23 - self.f123) * DNA23 + self.f123 * self.DNAtot) / self.gammamRNA

        # Set all variables of the steady state.
        self.ss = {'DNA0': DNA0, 'DNA1': DNA1, 'DNA2': DNA2, 'DNA3': DNA3, 'DNA12': DNA12, 'DNA13': DNA13,
                   'DNA23': DNA23, 'DNA123': DNA123, 'mRNA': mRNA, 'm': m, 'd': d}

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
        denom = (self.kb1 * self.kb2 * self.kb3 / (
            self.ku1 * self.ku2 * self.ku3) * self.bcoop12 * self.bcoop23 * self.bcoop13 * self.bcoop123 / (
                     self.ucoop12 * self.ucoop23 * self.ucoop13 * self.ucoop123) * d ** 3 + (self.kb1 * self.kb2 / (
            self.ku1 * self.ku2) * self.bcoop12 / self.ucoop12 + self.kb1 * self.kb3 / (
                                                                                                 self.ku1 * self.ku3) * self.bcoop13 / self.ucoop13 + self.kb2 * self.kb3 / (
                                                                                                 self.ku2 * self.ku3) * self.bcoop23 / self.ucoop23) * d ** 2 + (
                     self.kb1 / self.ku1 + self.kb2 / self.ku2 + self.kb3 / self.ku3) * d + 1)
        DNA0 = self.DNAtot / denom
        DNA1 = self.DNAtot * self.kb1 / self.ku1 * d / denom
        DNA2 = self.DNAtot * self.kb2 / self.ku2 * d / denom
        DNA3 = self.DNAtot * self.kb3 / self.ku3 * d / denom
        DNA12 = self.DNAtot * self.kb1 * self.kb2 / (
            self.ku1 * self.ku2) * self.bcoop12 / self.ucoop12 * d ** 2 / denom
        DNA13 = self.DNAtot * self.kb1 * self.kb3 / (
            self.ku1 * self.ku3) * self.bcoop13 / self.ucoop13 * d ** 2 / denom
        DNA23 = self.DNAtot * self.kb2 * self.kb3 / (
            self.ku2 * self.ku3) * self.bcoop23 / self.ucoop23 * d ** 2 / denom
        DNA123 = self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23
        qss = {'DNA0': DNA0, 'DNA1': DNA1, 'DNA2': DNA2, 'DNA3': DNA3, 'DNA12': DNA12, 'DNA13': DNA13, 'DNA23': DNA23,
              'DNA123': DNA123, 'd': d, 'm': m}

        return qss

    def jacobian(self, var):
        """ Return the jacobian of the system for variables = arg. """

        # Get the variables.
        DNA0, DNA1, DNA2, DNA3, DNA12, DNA23, DNA13, mRNA, m, d = var

        # Calculate some auxiliary parameters.
        c12DNA123 = self.ucoop123 * self.ucoop13 * self.ucoop23 * self.ku3
        c13DNA123 = self.ucoop123 * self.ucoop12 * self.ucoop23 * self.ku2
        c23DNA123 = self.ucoop123 * self.ucoop12 * self.ucoop13 * self.ku1
        cmDNA123 = self.ucoop123 * (
            self.ku1 * self.ucoop12 * self.ucoop13 + self.ku2 * self.ucoop12 * self.ucoop23 + self.ku3 * self.ucoop13 * self.ucoop23)

        # Define the Jacobian.
        J =  np.array(
            [[-(self.kb1 + self.kb2 + self.kb3) * d, self.ku1, self.ku2, self.ku3, 0, 0, 0, 0, 0,
              -(self.kb1 + self.kb2 + self.kb3) * DNA0],
              [self.kb1 * d, -(self.ku1 + (self.bcoop12 * self.kb2 + self.bcoop13 * self.kb3) * d), 0, 0,
               self.ucoop12 * self.ku2, self.ucoop13 * self.ku3, 0, 0, 0,
               self.kb1 * DNA0 - (self.bcoop12 * self.kb2 + self.bcoop13 * self.kb3) * DNA1],
              [self.kb2 * d, 0, -(self.ku2 + (self.bcoop12 * self.kb1 + self.bcoop23 * self.kb3) * d), 0,
               self.ucoop12 * self.ku1, 0, self.ucoop23 * self.ku3, 0, 0,
               self.kb2 * DNA0 - (self.bcoop12 * self.kb1 + self.bcoop23 * self.kb3) * DNA2],
              [self.kb3 * d, 0, 0, -(self.ku3 + (self.bcoop13 * self.kb1 + self.bcoop23 * self.kb2) * d),
               0, self.ucoop13 * self.ku1, self.ucoop23 * self.ku2, 0, 0,
               self.kb3 * DNA0 - (self.bcoop13 * self.kb1 + self.bcoop23 * self.kb2) * DNA3],
              [-c12DNA123, self.bcoop12 * d * self.kb2 - c12DNA123, self.bcoop12 * d * self.kb1 - c12DNA123,
               -c12DNA123, -(self.bcoop123 * self.bcoop13 * self.bcoop23 * self.kb3 * d + self.ucoop12 * (
                  self.ku1 + self.ku2)) - c12DNA123, -c12DNA123, -c12DNA123, 0, 0,
               self.bcoop12 * (self.kb1 * DNA2 + self.kb2 * DNA1) - (
                   self.bcoop123 * self.bcoop13 * self.bcoop23 * self.kb3) * DNA12],
              [-c13DNA123, self.bcoop13 * d * self.kb3 - c13DNA123, -c13DNA123,
               self.bcoop13 * d * self.kb1 - c13DNA123, -c13DNA123, -(
                  self.bcoop123 * self.bcoop12 * self.bcoop23 * self.kb2 * d + self.ucoop13 * (
                      self.ku1 + self.ku3)) - c13DNA123, -c13DNA123, 0, 0,
               self.bcoop13 * (self.kb1 * DNA3 + self.kb3 * DNA1) - (
                   self.bcoop123 * self.bcoop12 * self.bcoop23 * self.kb2) * DNA13],
              [-c23DNA123, -c23DNA123, self.bcoop23 * d * self.kb3 - c23DNA123,
               self.bcoop23 * d * self.kb2 - c23DNA123, -c23DNA123, -c23DNA123, -(
                  self.bcoop123 * self.bcoop12 * self.bcoop13 * self.kb1 * d + self.ucoop23 * (
                      self.ku2 + self.ku3)) - c23DNA123, 0, 0, self.bcoop23 * (self.kb3 * DNA2 + self.kb2 * DNA3) - (
                   self.bcoop123 * self.bcoop12 * self.bcoop13 * self.kb1) * DNA23],
              [self.phi0 * (1 - self.f123), self.phi0 * (self.f1 - self.f123), self.phi0 * (self.f2 - self.f123),
               self.phi0 * (self.f3 - self.f123), self.phi0 * (self.f12 - self.f123), self.phi0 * (self.f13 - self.f123),
               self.phi0 * (self.f23 - self.f123), -self.gammamRNA, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, self.beta, -4 * self.alphaass * m - self.gammam, 2 * self.alphadiss],
              [-d * (self.kb1 + self.kb2 + self.kb3) - cmDNA123,
               -d * (self.kb2 * self.bcoop12 + self.kb3 * self.bcoop13) + self.ku1 - cmDNA123,
               -d * (self.kb1 * self.bcoop12 + self.kb3 * self.bcoop23) + self.ku2 - cmDNA123,
               -d * (self.kb1 * self.bcoop13 + self.kb2 * self.bcoop23) + self.ku3 - cmDNA123,
               -d * self.kb3 * self.bcoop123 * self.bcoop13 * self.bcoop23 + (
                   self.ku1 + self.ku2) * self.ucoop12 - cmDNA123,
               -d * self.kb2 * self.bcoop123 * self.bcoop12 * self.bcoop23 + (
                   self.ku1 + self.ku3) * self.ucoop13 - cmDNA123,
               -d * self.kb1 * self.bcoop123 * self.bcoop12 * self.bcoop13 + (
                   self.ku2 + self.ku3) * self.ucoop23 - cmDNA123, 0, 2 * self.alphaass * m,
               -self.alphadiss - self.gammad - ((self.kb1 + self.kb2 + self.kb3) * DNA0 + (
                   self.kb2 * self.bcoop12 + self.kb3 * self.bcoop13) * DNA1 + (
                                                    self.kb1 * self.bcoop12 + self.kb3 * self.bcoop23) * DNA2
                                                + (self.kb1 * self.bcoop13 + self.kb2 * self.bcoop23) * DNA3
                                                + self.kb3 * self.bcoop123 * self.bcoop13 * self.bcoop23 * DNA12
                                                + self.kb2 * self.bcoop123 * self.bcoop12 * self.bcoop23 * DNA13
                                                + self.kb1 * self.bcoop123 * self.bcoop12 * self.bcoop13 * DNA23)]])

        return J

    # --------------------------------------------------------------------------------------------------------
    # AUXILIARY FUNCTIONS
    # --------------------------------------------------------------------------------------------------------

    def dimensionless_parameters(self, ts, qs):
        """ Return parameters without dimensions using time scale ts and quantity scale qs."""

        par = {  # nondimensionalized
            'beta': self.beta * ts, 'gammam': self.gammam * ts, 'gammamRNA': self.gammamRNA * ts,
            'gammad': self.gammad * ts, 'alphaass': self.alphaass * ts * qs, 'alphadiss': self.alphadiss * ts,
            'taum': self.taum / ts, 'taumRNA': self.taumRNA / ts, 'DNAtot': self.DNAtot / qs, 'phi0': self.phi0 * ts,
            'fa': self.fa, 'fr': self.fr, 'kb1': self.kb1 * qs * ts, 'ku1': self.ku1 * ts,
            'kb2': self.kb2 * qs * ts, 'ku2': self.ku2 * ts, 'kb3': self.kb3 * qs * ts, 'ku3': self.ku3 * ts,
            'bcoop12': self.bcoop12, 'bcoop23': self.bcoop23, 'bcoop13': self.bcoop13, 'bcoop123': self.bcoop123,
            'ucoop12': self.ucoop12, 'ucoop23': self.ucoop23, 'ucoop13': self.ucoop13, 'ucoop123': self.ucoop123,
            'volume': self.vol, 'ts': ts, 'qs': qs}
        return par

    def SsLrpB_compatible(self):
        """
        Checks whether the system is Ss-LrpB compatible, i.e. the parameters are compatible with the measured ones,
        the non-monotonicity type is 1 and the maximum response is higher than 2 times the basal rate
        and the minimum lower than 0.5 times the basal rate.
        """

        xx = np.linspace(0, 1000, 1000)

        y = self.transcription_rate(xx, polymer=Polymer.DIMER)

        ct = 1e-3 / 2.4

        right_parameters = [np.abs((self.kb1 / self.ku1 - 73.5 * ct) / 73.5 * ct) < 1e-5,# and
                            np.abs((self.kb2 / self.ku2 - 0.7 * ct) / 0.7 * ct) < 1e-5,# and
                            np.abs((self.kb3 / self.ku3 - 49.1 * ct) / 49.1 * ct) < 1e-5,# and
                            np.abs((self.bcoop12 / self.ucoop12 - 4.1) / 4.1) < 1e-2,# and
                            np.abs((self.bcoop23 / self.ucoop23 - 2.1) / 2.1) < 1e-2,# and
                            np.abs((self.bcoop13 / self.ucoop13 - 7.9) / 7.9) < 1e-2,# and
                            np.abs((self.bcoop123 / self.ucoop123 - 3.1) / 3.1) < 1e-2]

        return right_parameters and self.monotonicity_type() == 1 and max(y) / y[0] > 2 and min(y) / y[0] < 0.5

