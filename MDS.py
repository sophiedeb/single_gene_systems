#!/usr/bin/env python

"""MDS.py: Class monomer dimer system (MDS)."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

from singlegenesystem import *

class MDS(SG):
    DNA_STATES = np.array(['DNA0', 'DNAm', 'DNAd'])
    ALL_VARS = np.array(['DNA0', 'DNAm', 'mRNA', 'm', 'd'])
    PARAMETERS = np.array(
        ['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm', 'fd', 'Km', 'Kd', 'kbm', 'kbd',
         'DNAtot', 'taum', 'taumRNA'])

    def __init__(self, params={}):
        self.ss = {}  # steady state

        self.set_parameters(params)

    # --------------------------------------------------------------------------------------------------------
    # PARAMETER FUNCTIONS
    # --------------------------------------------------------------------------------------------------------

    def set_parameters(self, params):
        """ Set parameters of the MDS.

        Only parameters of parameterlist are accepted, other parameters are initialized by None.
        Automatic conversion from binding constants to unbinding rates and from cooperativity to unbinding cooperativity.
        The steady state is cleared, because different parameters will lead to a different steady state.
        """

        for param in self.PARAMETERS:
            if param in params:
                setattr(self, param, params[param])
            else:
                setattr(self, param, None)

        for i in ['m', 'd']:
            if 'K%c' % i in params:
                setattr(self, 'ku%c' % i, params['kb%c' % i] / params['K%c' % i])
            elif 'ku%c' % i in params:
                setattr(self, 'K%c' % i, params['kb%c' % i] / params['ku%c' % i])
                setattr(self, 'ku%c' % i, params['ku%c' % i])

        self.ss.clear()

    def __copy__(self):
        """ Makes a copy."""

        newone = type(self)()
        newone.set_parameters(self.get_parameters())

        return newone

    # --------------------------------------------------------------------------------------------------------
    # EQUATIONS AND TIME SERIES
    # --------------------------------------------------------------------------------------------------------

    def eqns(self, var, t, fastdimer=False):
        """ Ordinary differential equations of system.
        """

        if not fastdimer:
            # Set variables
            DNA0, DNAm, mRNA, m, d = var

            # Calculate all differential equations.
            dDNA0 = -(self.kbm * m + self.kbd * d) * DNA0 + self.kum * DNAm + self.kud * (self.DNAtot - DNA0 - DNAm)
            dDNAm = self.kbm * m * DNA0 - self.kum * DNAm
            dmRNA = self.phi0 * (DNA0 + self.fm * DNAm + self.fd * (self.DNAtot - DNA0 - DNAm)) - self.gammamRNA * mRNA
            dm = -self.kbm * m * DNA0 + self.kum * DNAm + self.beta * mRNA - 2 * self.alphaass * pow(m,
                                                                                                     2) + 2 * self.alphadiss * d - self.gammam * m
            dd = -self.kbd * d * DNA0 + self.kud * (self.DNAtot - DNA0 - DNAm) + self.alphaass * pow(m,
                                                                                                     2) - self.alphadiss * d - self.gammad * d
            return [dDNA0, dDNAm, dmRNA, dm, dd]
        else:
            # Set variables
            mRNA, m = var

            # Calculate all differential equations.
            d = self.alphaass * m ** 2 / (self.alphadiss + self.gammad)
            DNA0 = (-self.Kd * d + self.DNAtot) / (1 + self.Km * m)
            DNAm = (-self.Kd * d + self.DNAtot) / (1 + 1 / (self.Km * m))
            dmRNA = self.phi0 * DNA0 + self.phi0 * self.fm * DNAm + self.phi0 * self.fd * (
                self.DNAtot - DNA0 - DNAm) - self.gammamRNA * mRNA
            dm = - self.kbm * m * DNA0 + self.kum * DNAm + self.beta * mRNA - 2 * self.alphaass * pow(m,
                                                                                                      2) + 2 * self.alphadiss * d - self.gammam * m
            return [dmRNA, dm]

    # delay differential equations in dictionary format for pydelay
    @staticmethod
    def eqnsdelay(fastdimer=False):
        """ Delay differential equations.

        Use of a dictionary of functions for package pydelay.
        """

        if not fastdimer:
            eq = {'DNA0': '-(kbm*m+kbd*d)*DNA0 + kum*DNAm + kud*(DNAtot-DNA0-DNAm)', 'DNAm': 'kbm*m*DNA0 - kum*DNAm',
                  'mRNA': 'phi0*DNA0 + phi0*fm*DNAm + phi0*fd*(DNAtot-DNA0-DNAm) -gammamRNA*mRNA',
                  'm': '- kbm*m*DNA0 + kum*DNAm + beta*mRNA - 2*alphaass*pow(m,2) + 2*alphadiss*d - gammam*m',
                  'd': '-kbd*d*DNA0 + kud*(DNAtot-DNA0-DNAm) + alphaass*pow(m,2) - alphadiss*d - gammad*d'}
        else:
            eq = {'m': 'beta*mRNA(t-taum) - 2*alphaass*pow(m,2)/(alphadiss/gammad+1) - gammam*m',
                  'mRNA': '(phi0*(1-fr))*DNAtot*(alphadiss + gammad)/(Kd*alphaass*pow(m,2)'
                          '+ (Km*m+1)*(alphadiss+gammad)) +(phi0*(fa-fr))*DNAtot*Km*m*(alphadiss + gammad)/'
                          '(Kd*alphaass*pow(m,2) + (Km*m+1)*(alphadiss + gammad))'
                          '+phi0*fr*DNAtot'
                          '-gammamRNA*mRNA'}
        return eq

    @staticmethod
    def reaction_matrix():
        """ Returns the reaction matrix.

        Columns are chemical processes. Rows are the change of every variable at given process.
        """

        # Columns are chemical processes. For order see function propensities.
        # Rows are the change of the variables after specific chemical process.
        # Rows in order : DNA0,
        reactions = np.matrix([[-1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0], [-1, 0, 1, 0, 0, 0, 1, -2, 2, -1, 0],
                               [0, -1, 0, 1, 0, 0, 0, 1, -1, 0, -1]])
        return reactions

    def propensities(self, var):
        """ Propensities for variables = arg."""

        # Get variables.
        DNA0, DNAm, mRNA, m, d = var
        DNAd = self.DNAtot - DNA0 - DNAm

        # Set propensities for all processes.
        prop = np.array([self.kbm * DNA0 * m,  # DNA0 + m --> DNAm (binding)
                         self.kbd * d * DNA0,  # DNA0 + d --> DNAd (binding)
                         self.kum * DNAm,  # DNAm --> DNA0 + m (unbinding)
                         self.kud * DNAd,  # DNAd --> DNA0 + d (unbinding)
                         self.phi0 * (DNA0 + self.fm * DNAm + self.fd * DNAd),  # --> mRNA (transcription)
                         self.gammamRNA * mRNA,  # mRNA --> 0 (mRNA degradation)
                         self.beta * mRNA,  # mRNA --> mRNA + m (translation)
                         self.alphaass * m * (m - 1),  # 2*m --> d (dimer association)
                         self.alphadiss * d,  # d --> 2*m (dimer dissociation)
                         self.gammam * m,  # m --> 0 (monomer degradation)
                         self.gammad * d])  # d --> 0 (dimer degradation)

        # Print message if any propensity is negative.
        if np.any(prop < 0):
            print(prop)

        return prop

    # --------------------------------------------------------------------------------------------------------
    # STEADY STATE, JACOBIAN & EIGENVALUES
    # --------------------------------------------------------------------------------------------------------

    def set_steady_state(self):
        """ Set the steady state of the system."""

        # Define all auxiliary functions.
        fd = lambda x: self.alphaass * x ** 2 / (self.alphadiss + self.gammad)
        fDNA0 = lambda x: self.DNAtot / (self.Kd * fd(x) + self.Km * x + 1)
        fDNAm = lambda x: self.Km * x * fDNA0(x)
        fmRNA = lambda x: self.phi0 * (
            (1 - self.fd) * fDNA0(x) + (self.fm - self.fd) * fDNAm(x) + self.fd * self.DNAtot) / self.gammamRNA

        # Function of the monomer concentration for dmRNA/dt = 0.
        func = lambda x: ((self.phi0 * (1 - self.fd)) * (self.DNAtot / (self.Kd * fd(x) + self.Km * x + 1)) + (
            self.phi0 * (self.fm - self.fd)) * (self.DNAtot * self.Km * x / (
            self.Kd * fd(x) + self.Km * x + 1)) + self.phi0 * self.fd * self.DNAtot - self.gammamRNA / self.beta * (
                              self.gammam * x) - self.gammamRNA / self.beta * 2 * (
                              self.alphaass * x ** 2 - self.alphadiss * fd(x)))

        # Find steady state(s) of the monomer.
        m = self._find_x_steady_state(func)

        # Calculate the steady state values for every variable from the monomer steady state.
        d = fd(m)
        DNA0 = fDNA0(m)
        DNAm = fDNAm(m)
        mRNA = fmRNA(m)

        # Set all variables of the steady state.
        self.ss = {'DNA0': DNA0, 'DNAm': DNAm, 'DNAd': self.DNAtot - DNA0 - DNAm, 'mRNA': mRNA, 'm': m, 'd': d}

        return self.ss

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
        DNA0 = self.DNAtot / (self.Kd * d + self.Km * m + 1)
        DNAm = self.Km * m * DNA0
        DNAd = self.Kd * d * DNA0
        qss = {'DNA0': DNA0, 'DNAm': DNAm, 'DNAd': DNAd, 'd': d, 'm': m}

        return qss


    def jacobian(self, var, fastdimer=False):
        """ Return the jacobian of the system for variables = arg. """

        # Get the variables.
        DNA0, DNAm, mRNA, m, d = var

        # Define the Jacobian.
        if not fastdimer:
            J = [
                [-(self.kbm * m + self.kbd * d + self.kud), self.kum - self.kud, 0, -self.kbm * DNA0, -self.kbd * DNA0],
                [self.kbm * m, -self.kum, 0, self.kbm * DNA0, 0],
                [self.phi0 * (1 - self.fd), self.phi0 * (self.fm - self.fd), -self.gammamRNA, 0, 0],
                [-self.kbm * m, self.kum, self.beta, -self.kbm * DNA0 - 4 * self.alphaass * m - self.gammam,
                 2 * self.alphadiss], [-self.kbd * d - self.kud, -self.kud, 0, 2 * self.alphaass * m,
                                       -(self.kbd * DNA0 + self.alphadiss + self.gammad)]]
        else:
            DNA0 = self.DNAtot * (self.alphadiss + self.gammad) / (
                self.Kd * self.alphaass * m ** 2 + (self.Km * m + 1) * (self.alphadiss + self.gammad))
            dDNA0 = -self.DNAtot * (self.alphadiss + self.gammad) * (
                2 * self.Kd * self.alphaass * m + self.Km * (self.alphadiss + self.gammad) / (
                    self.Kd * self.alphaass * m ** 2 + (self.Km * m + 1) * (self.alphadiss + self.gammad)))
            J = [[-4 * self.alphaass * m / (self.alphadiss + self.gammad) - self.gammam, self.beta], [
                (self.phi0 * (1 - self.fd)) * dDNA0 + (self.phi0 * (self.fm - self.fd)) * (
                    self.Km * m * dDNA0 + self.Km * DNA0), -self.gammamRNA]]

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
            'fm': self.fm, 'fd': self.fd, 'Km': self.Km * qs, 'Kd': self.Kd * qs, 'kbm': self.kbm * qs * ts,
            'kum': self.kum * ts, 'kbd': self.kbd * qs * ts, 'kud': self.kud * ts, 'volume': self.vol, 'ts': ts,
            'qs': qs}

        return par
