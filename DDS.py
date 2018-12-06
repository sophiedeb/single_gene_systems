from singlegenesystem import *

class DDS(SG):
    def __init__(self, params={}):
        # self.par = params;
        if(len(params)==0):
            self.beta = 0
            self.gammam = 0
            self.gammamRNA = 0
            self.gammad = 0
            self.alphaass = 0
            self.alphadiss = 0
            self.taum = 0
            self.taumRNA = 0
            self.DNAtot = 0
            self.phi0 = 0
            self.f1 = 0
            self.f2 = 0
            self.f12 = 0
            self.kb1 = 0
            self.ku1 = 0
            self.Kd1 = 0
            self.kb2 = 0
            self.ku2 = 0
            self.Kd2 = 0
            self.bcoop = 0
            self.ucoop = 0
            self.vol = 0
        else:
            self.beta = params['beta']
            self.gammam = params['gammam'];
            self.gammamRNA = params['gammamRNA'];
            self.gammad = params['gammad'];
            self.alphaass = params['alphaass'];
            if('alphadiss' in params):
                self.alphadiss = params['alphadiss'];
                self.Kdim = self.alphadiss/self.alphaass;
            elif('Kdim' in params):
                self.Kdim = params['Kdim'];
                self.alphadiss = self.alphaass*self.Kdim;
            self.taum = params['taum'];
            self.taumRNA = params['taumRNA'];
            self.DNAtot = params['DNAtot'];
            self.phi0 = params['phi0'];
            self.f1 = params['f1'];
            self.f2 = params['f2'];
            self.f12 = params['f12'];
            self.kb1 = params['kb1'];
            self.ku1 = params['ku1'];
            self.Kd1 = self.kb1 / self.ku1;
            self.kb2 = params['kb2'];
            self.ku2 = params['ku2'];
            self.Kd2 = self.kb2 / self.ku2;
            self.bcoop = params['bcoop'];
            self.ucoop = params['ucoop'];
            self.vol = params['vol'];
        self.ss = {}
        self.DNAstates = np.array(['DNA0', 'DNA1', 'DNA2', 'DNA12']);
        self.allvars = np.array(['DNA0', 'DNA1', 'DNA2', 'mRNA', 'm', 'd']);
        self.nameparameters = np.array(['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', #TODO fa, fr -> f1, f2, f12
            'f2', 'f12', 'kb1', 'ku1', 'kb2', 'ku2', 'bcoop', 'ucoop'])

    def outRange(self):
        physrange = physiologicalRange();
        line = "";
        par = [self.beta, self.gammam, self.gammamRNA, self.gammad, self.beta, self.alphaass, self.alphadiss, self.phi0,
               self.f1, self.f2, self.f12, self.kb1, self.ku1, self.kb2, self.ku2, self.bcoop, self.ucoop];
        name = ['beta', 'gammam', 'gammamRNA', 'gammad', 'beta', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f12', 'kb1',
                   'ku1', 'kb2', 'ku2', 'bcoop', 'ucoop']
        parname = ['beta', 'gammam', 'gammamRNA', 'gammad', 'beta', 'alphaass', 'alphadiss', 'phi0', 'fa', 'fa', 'fr', 'kb',
                   'ku', 'kb', 'ku', 'bcoop', 'ucoop']
        for i in range(len(par)):
            if(par[i] < physrange[parname[i]][0]):
                line += "%s low; "%name[i];
            elif(par[i] == physrange[parname[i]][0]):
                line += "%s low bound; " % name[i];
            elif (par[i] > physrange[parname[i]][1]):
                line += "%s high; " % name[i];
            elif (par[i] == physrange[parname[i]][1]):
                line += "%s high bound; " % name[i];
        return line;

    def nullClineM(self, x, monomer=True):
        qss = self.quasisteadystate(x, monomer);
        # values of mRNA for which dm/dt are 0 :
        m0 = (2 * self.alphaass * qss['m'] ** 2 - 2 * self.alphadiss * qss['d'] + self.gammam * qss['m']) / self.beta;
        return m0;

    def nullClineMRNA(self, x, monomer=True):
        qss = self.quasisteadystate(x, monomer);
        # values of mRNA for which dmRNA/dt are 0 :
        mRNA0 = self.phi0 * (
            qss['DNA0'] + self.f1 * qss['DNA1'] + self.f2 * qss['DNA2'] + self.f12 * (self.DNAtot - qss['DNA0'] - qss['DNA1'] - qss['DNA2'])) / self.gammamRNA;
        return mRNA0;

    def parameterline(self):
        string = "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E \n" % (
            self.beta, self.gammam, self.gammamRNA, self.gammad, self.alphaass, self.alphadiss, self.phi0, self.f1, #TODO fa, fr -> f1, f2, f12
            self.f2, self.f12, self.kb1, self.ku1, self.kb2, self.ku2, self.bcoop, self.ucoop);
        return string

    def copy(self, parent):
        self.beta = parent.beta;
        self.gammam = parent.gammam;
        self.gammamRNA = parent.gammamRNA;
        self.gammad = parent.gammad;
        self.alphaass = parent.alphaass;
        self.Kdim = parent.Kdim;
        self.alphadiss = self.Kdim * self.alphaass;
        self.taum = parent.taum;
        self.taumRNA = parent.taumRNA;
        self.DNAtot = parent.DNAtot;
        self.phi0 = parent.phi0;
        self.f1 = parent.f1
        self.f2 = parent.f2
        self.f12 = parent.f12
        self.kb1 = parent.kb1
        self.kb2 = parent.kb2
        self.ku1 = parent.ku1
        self.ku2 = parent.ku2
        self.Kd1 = self.kb1 / self.ku1;
        self.Kd2 = self.kb2 / self.ku2;
        self.bcoop = parent.bcoop
        self.ucoop = parent.ucoop

    def smallmutant(self, parent, constraints={}, SsLrpB=False):
        idx = int(rd.uniform(0, 16))
        mutation = np.ones(16);
        mutation[idx] = 10**(rd.uniform(np.log10(0.5), np.log10(2)));
        if(len(constraints)==0):
            print("no constraints")
            self.beta = parent.beta * mutation[0];
            self.gammam = parent.gammam * mutation[1];
            self.gammamRNA = parent.gammamRNA * mutation[2];
            self.gammad = parent.gammad * mutation[3];
            self.alphaass = parent.alphaass * mutation[4];
            self.Kdim = parent.Kdim * mutation[5];
            self.alphadiss = self.Kdim * self.alphaass;  # parent.alphadiss * mutation[5];
            self.taum = parent.taum;
            self.taumRNA = parent.taumRNA;
            self.DNAtot = parent.DNAtot;
            self.phi0 = parent.phi0 * mutation[6];
            self.fa = parent.fa * mutation[7];
            self.fr = parent.fr * mutation[8];
            self.kb1 = parent.kb1 * mutation[9];
            self.ku1 = parent.ku1 * mutation[10];
            self.kb2 = parent.kb2 * mutation[11];
            self.ku2 = parent.ku2 * mutation[12];
            self.Kd1 = self.kb1 / self.ku1;
            self.Kd2 = self.kb2 / self.ku2;
            self.bcoop = parent.bcoop * mutation[13];
            self.ucoop = parent.ucoop * mutation[14];
            self.vol = parent.vol;
        else:
            self.beta = max(min(parent.beta * mutation[0], constraints['beta'][1]), constraints['beta'][0]);
            self.gammam = max(min(parent.gammam * mutation[1], constraints['gammam'][1]), constraints['gammam'][0]);
            self.gammamRNA = max(min(parent.gammamRNA * mutation[2], constraints['gammamRNA'][1]), constraints['gammamRNA'][0]);
            self.gammad = max(min(parent.gammad * mutation[3], constraints['gammad'][1]), constraints['gammad'][0]);
            self.alphaass = max(min(parent.alphaass * mutation[4], constraints['alphaass'][1]), constraints['alphaass'][0]);
            self.Kdim = max(min(parent.Kdim * mutation[5], constraints['Kdim'][1]), constraints['Kdim'][0]);
            self.alphadiss = self.Kdim * self.alphaass;  # parent.alphadiss * mutation[5];
            #self.alphadiss = max(min(parent.alphadiss * mutation[5], constraints['alphadiss'][1]), constraints['alphadiss'][0]);
            self.taum = parent.taum;
            self.taumRNA = parent.taumRNA;
            self.DNAtot = parent.DNAtot;
            self.phi0 = max(min(parent.phi0 * mutation[6], constraints['phi0'][1]), constraints['phi0'][0]);
            self.f1 = max(min(parent.f1 * mutation[7], constraints['f'][1]), constraints['f'][0]);#todo
            self.f2 = max(min(parent.f2 * mutation[8], constraints['f'][1]), constraints['f'][0]);#todo
            self.f12 = max(min(parent.f2 * mutation[8], constraints['f'][1]), constraints['f'][0]);#todo #
            #self.f12 = max(min(parent.f12 * mutation[9], max(self.f1, self.f2)), constraints['f'][0]);#todo
            self.kb1 = max(min(parent.kb1 * mutation[10], constraints['kb'][1]), constraints['kb'][0]);
            self.ku1 = max(min(parent.ku1 * mutation[11], constraints['ku'][1]), constraints['ku'][0]);
            self.kb2 = max(min(parent.kb2 * mutation[12], constraints['kb'][1]), constraints['kb'][0]);
            self.ku2 = max(min(parent.ku2 * mutation[13], constraints['ku'][1]), constraints['ku'][0]);
            self.Kd1 = self.kb1 / self.ku1;
            self.Kd2 = self.kb2 / self.ku2;
            self.bcoop = max(min(parent.bcoop * mutation[14], constraints['bcoop'][1]), constraints['bcoop'][0]);
            self.ucoop = max(min(parent.ucoop * mutation[15], constraints['ucoop'][1]), constraints['ucoop'][0]);
            self.vol = parent.vol;
        self.ss.clear();
        return

    def largemutant(self, parent, constraints={}, SsLrpB=False):
        idx = int(rd.uniform(0, 16))
        self.copy(parent)
        if(idx==0):
            self.beta = np.exp(rd.uniform(np.log(constraints['beta'][0]), np.log(constraints['beta'][1])))
        elif (idx == 1):
            self.gammam = np.exp(rd.uniform(np.log(constraints['gammam'][0]), np.log(constraints['gammam'][1])))
        elif (idx == 2):
            self.gammamRNA = np.exp(rd.uniform(np.log(constraints['gammamRNA'][0]), np.log(constraints['gammamRNA'][1])))
        elif (idx == 3):
            self.gammad = np.exp(rd.uniform(np.log(constraints['gammad'][0]), np.log(constraints['gammad'][1])))
        elif (idx == 4):
            self.alphaass = np.exp(rd.uniform(np.log(constraints['alphaass'][0]), np.log(constraints['alphaass'][1])))
        elif (idx == 5):
            self.Kdim = np.exp(rd.uniform(np.log(constraints['Kdim'][0]), np.log(constraints['Kdim'][1])))
            self.alphadiss = self.Kdim * self.alphaass;
        elif (idx == 6):
            self.phi0 = np.exp(rd.uniform(np.log(constraints['phi0'][0]), np.log(constraints['phi0'][1])))
        elif (idx == 7):
            self.f1 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo
            #self.f12 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(max(self.f1,self.f2))))#todo
        elif (idx == 8):
            self.f2 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1])))#todo
            #self.f12 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(max(self.f1,self.f2))))#todo
        elif (idx == 9):
            self.f12 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1])))#todo
            #self.f12 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(max(self.f1,self.f2))))#todo
        elif (idx == 10):
            self.kb1 = np.exp(rd.uniform(np.log(constraints['kb'][0]), np.log(constraints['kb'][1])))
        elif (idx == 11):
            self.kb2 = np.exp(rd.uniform(np.log(constraints['kb'][0]), np.log(constraints['kb'][1])))
        elif (idx == 12):
            self.ku1 = np.exp(rd.uniform(np.log(constraints['ku'][0]), np.log(constraints['ku'][1])))
        elif (idx == 13):
            self.ku2 = np.exp(rd.uniform(np.log(constraints['ku'][0]), np.log(constraints['ku'][1])))
        elif (idx == 14):
            self.bcoop = np.exp(rd.uniform(np.log(constraints['bcoop'][0]), np.log(constraints['bcoop'][1])))
        elif (idx == 15):
            self.ucoop = np.exp(rd.uniform(np.log(constraints['ucoop'][0]), np.log(constraints['ucoop'][1])))
        self.Kd1 = self.kb1 / self.ku1;
        self.Kd2 = self.kb2 / self.ku2;
        self.ss.clear();
        return

    def sscurve(self, ax, e=5000, c='k'):
        if (not len(self.ss) > 0):
            self.setsteadystate();

        fd = lambda x: self.alphaass * x ** 2 / (self.alphadiss + self.gammad);
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

    def setsteadystate(self):
        fd = lambda x: self.alphaass * x ** 2 / (self.alphadiss + self.gammad);
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

        if (False):
            x = np.logspace(-7, 3, 2000)
            #x = np.linspace(0, 100, 2000)
            plt.figure('%.3E %.3E %.3E'%(self.f1, self.f2, self.f12))
            plt.plot(x, func(x))
            plt.xscale('log')
            plt.ylim([-100, 100])
            # plt.plot(m, 0, '*')
            plt.show()
            print(np.sign(func(1e25))==np.sign(func(1e-25)), np.sign(func(1e25)), np.sign(func(1e-25)))
        if (func(1e25) > 0):
            m = np.array([np.inf])
        elif (func(1e-25) < 0):
            m = np.array([0])
        elif (np.sign(func(1e25))==np.sign(func(1e-25))):
            m = np.nan
        else:
            N=5000
            x = np.logspace(-25,25,N);
            zeroes = np.where(np.diff(np.sign(func(x))))[0];
            if(len(zeroes) > 1): # multistability
                zeroes = np.append(np.append([0], zeroes), [N-1])
                m = np.array([brentq(func, (x[zeroes[i]]+x[zeroes[i-1]+1])/2,
                                     (x[zeroes[i+1]]+x[zeroes[i]+1])/2,
                                xtol=1e-10, rtol=1e-5, maxiter=200) for i in range(1,len(zeroes)-1)])
            else:
                m = np.array([brentq(func, 1e-25, 1e25, xtol=1e-10, rtol=1e-5, maxiter=200)])

        d = fd(m)
        DNA0 = f0(m)
        DNA1 = f1(m)
        DNA2 = f2(m)
        mRNA = self.phi0 * ((1 - self.f12)*DNA0 + (self.f1 - self.f12) * DNA1 + (self.f2 - self.f12) * DNA2 + self.f12 * self.DNAtot) / self.gammamRNA
        self.ss = {'DNA0': DNA0, 'DNA1': DNA1, 'DNA2': DNA2, 'DNA12': self.DNAtot - DNA0 - DNA1 - DNA2, 'mRNA': mRNA,
                   'm': m, 'd': d}

    def setparameters(self, params):
        self.beta = params['beta']
        self.gammam = params['gammam'];
        self.gammamRNA = params['gammamRNA'];
        self.gammad = params['gammad'];
        self.alphaass = params['alphaass'];
        if ('alphadiss' in params):
            self.alphadiss = params['alphadiss'];
            self.Kdim = self.alphadiss / self.alphaass;
        elif ('Kdim' in params):
            self.Kdim = params['Kdim'];
            self.alphadiss = self.alphaass * self.Kdim;
        self.taum = params['taum'];
        self.taumRNA = params['taumRNA'];
        self.DNAtot = params['DNAtot'];
        self.phi0 = params['phi0'];
        self.f1 = params['f1'];
        self.f2 = params['f2']
        self.f12 = params['f12']
        self.kb1 = params['kb1'];
        if ('ku1' in params):
            self.ku1 = params['ku1'];
            self.Kd1 = self.kb1 / self.ku1;
        else:
            self.Kd1 = params['Kd1'];
            self.ku1 = self.kb1 / self.Kd1;
        self.kb2 = params['kb2'];
        if ('ku2' in params):
            self.ku2 = params['ku2'];
            self.Kd2 = self.kb2 / self.ku2;
        else:
            self.Kd2 = params['Kd2'];
            self.ku2 = self.kb2 / self.Kd2;
        self.bcoop = params['bcoop'];
        self.ucoop = params['ucoop'];
        self.vol = params['vol'];
        self.ss.clear()

    def steadystate(self):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        return self.ss

    def eqns(self, var, t=0.0, ts=1.0, qs=1.0):
        DNA0, DNA1, DNA2, mRNA, m, d = var;
        # ts, qs = args;

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
        return [dDNA0, dDNA1, dDNA2, dmRNA, dm, dd];

    def eqnsreduc(self, var, t=0.0, ts=1.0, qs=1.0):
        mRNA, m = var;
        # ts, qs = args;

        qss = self.quasisteadystate(m)
        DNA0 = qss['DNA0']
        DNA1 = qss['DNA1']
        DNA2 = qss['DNA2']
        d = qss['d']

        dmRNA = self.phi0 * ((1 - self.f12)*DNA0 + (self.f1 - self.f12) * DNA1 + (self.f2 - self.f12) * DNA2 + self.f12 * self.DNAtot) - self.gammamRNA * mRNA
        dm = self.beta * mRNA - 2 * self.alphaass * pow(m, 2) + 2 * self.alphadiss * d - self.gammam * m
        return [dmRNA, dm];

    def reactionmatrix(self):
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
        DNA0, DNA1, DNA2, mRNA, m, d = var;

        prop = np.array(
                  [self.kb1*DNA0*d,
                  self.ku1*DNA1,
                  self.kb2*DNA0*d,
                  self.ku2*DNA2,
                  self.kb2*self.bcoop*DNA1*d,
                  self.ku2*self.ucoop*(self.DNAtot-DNA0-DNA1-DNA2),
                  self.kb1*self.bcoop*DNA2*d,
                  self.ku1*self.ucoop*(self.DNAtot-DNA0-DNA1-DNA2),
                  self.phi0*((1 - self.f12)*DNA0 + (self.f1 - self.f12)*DNA1 + (self.f2 - self.f12)*DNA2 + self.f12*self.DNAtot),
                  self.beta*mRNA,
                  self.alphaass*m*(m-1),
                  self.alphadiss*d,
                  self.gammamRNA*mRNA,
                  self.gammam*m,
                  self.gammad*d])
        if(np.any(prop<0)):
            print(prop)
        return prop

    def histfunc(self, qs=1.0, fastdimer=False):
        if (fastdimer):
            hf = {'m': lambda t: self.ss['m'] / qs * 1.1, 'mRNA': lambda t: self.ss['mRNA'] / qs * 0.9}
        else:
            if (not len(self.ss) > 0):
                self.setsteadystate();
            hf = {'m': lambda t: self.ss['m'][0] / qs , 'd': lambda t: self.ss['d'][0] / qs * 3,
                  'mRNA': lambda t: self.ss['mRNA'][0] / qs , 'DNA0': lambda t: self.ss['DNA0'][0] / qs,
                  'DNA1': lambda t: self.ss['DNA1'][0] / qs, 'DNA2': lambda t: self.ss['DNA2'][0] / qs}
        return hf

    def quasisteadystate(self, x, monomer=True):
        if (monomer):
            m = x;
            d = self.alphaass * m ** 2 / (self.alphadiss + self.gammad)
        else:
            d = x;
            m = np.sqrt((self.alphadiss + self.gammad) * d / self.alphaass);

        fd = lambda x: self.alphaass * x ** 2 / (self.alphadiss + self.gammad);
        f0 = lambda x: self.DNAtot / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)
        f1 = lambda x: self.DNAtot * fd(x) * self.Kd1 / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)
        f2 = lambda x: self.DNAtot * fd(x) * self.Kd2 / (
            self.bcoop / self.ucoop * fd(x) ** 2 * self.Kd1 * self.Kd2 + fd(x) * (self.Kd1 + self.Kd2) + 1)

        DNA0 = f0(m)
        DNA1 = f1(m)
        DNA2 = f2(m)
        ss = {'DNA0': DNA0, 'DNA1': DNA1, 'DNA2': DNA2, 'DNA12': self.DNAtot - DNA0 - DNA1 - DNA2, 'd': d, 'm': m}
        return ss

    def transfunc(self):
        a = (self.Kd1*self.Kd2*self.bcoop/self.ucoop*(self.Kd1*(self.f12 - self.f1)+ self.Kd2*(self.f12 - self.f2)))
        b = (2*self.Kd1*self.Kd2*self.bcoop/self.ucoop*(self.f12 - 1))
        c = self.Kd1*(self.f1 - 1) + self.Kd2*(self.f2 - 1)
        D = b**2 - 4*a*c
        x1 = (-b-np.sqrt(D))/(2*a) if D>0 else np.nan
        x2 = (-b+np.sqrt(D))/(2*a) if D>0 else np.nan
        return D, x1, x2

    def transcriptionrate(self, m, monomer=True):
        ssDNA = self.quasisteadystate(m, monomer)
        return self.phi0 * (ssDNA['DNA0'] + self.f1 * ssDNA['DNA1'] + self.f2 * ssDNA['DNA2'] + self.f12 * ssDNA['DNA12'])

    def translationrate(self, m, monomer=True):
        return self.beta*self.transcriptionrate(m, monomer)/self.gammamRNA

    def paramsnd(self, ts, qs):
        par = {  # nondimensionalized
            'beta': self.beta * ts, 'gammam': self.gammam * ts, 'gammamRNA': self.gammamRNA * ts,
            'gammad': self.gammad * ts, 'alphaass': self.alphaass * ts * qs, 'alphadiss': self.alphadiss * ts,
            'taum': self.taum / ts, 'taumRNA': self.taumRNA / ts, 'DNAtot': self.DNAtot / qs, 'phi0': self.phi0 * ts,
            'fa': self.fa, 'fr': self.fr, 'Km': self.Km * qs, 'Kd': self.Kd * qs, 'kbm': self.kbm * qs * ts,
            'kum': self.kum * ts, 'kb': self.kb * qs * ts, 'ku': self.ku * ts, 'volume': self.vol, 'ts': ts, 'qs': qs}
        return par

    def timelapsereduced(self, t_f, dt, t_i=0, hf={}):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        ts = 1.0;
        qs = 1.0
        if (self.taum == 0 and self.taumRNA == 0):
            if (len(hf) == 0):  # no initial conditions given, start from perturbation of steady state
                if (not len(self.ss) > 0):
                    self.setsteadystate();
                initcond = [self.ss['mRNA'][0]*2, self.ss['m'][0] * 2]
            else:
                initcond = [hf[i](0) for i in ['mRNA', 'm']];
            t = np.arange(t_i, t_f, dt)
            ode = ode23(self, odeint(self.eqnsreduc, initcond, t, args=(ts, qs)), t, []);
            return ode
        return np.nan

    def timelapse(self, t_f, dt, t_i=0, case=0, hf={}):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        ts = 1.0;
        qs = 1.0
        if (self.taum == 0 and self.taumRNA == 0):
            if (len(hf) == 0):  # no initial conditions given, start from perturbation of steady state
                if (not len(self.ss) > 0):
                    self.setsteadystate();
                initcond = [self.ss['DNA0'][0], self.ss['DNA1'][0], self.ss['DNA2'][0], self.ss['mRNA'][0] * 1.05,
                            self.ss['m'][0] * 1.05, 3 * self.ss['d'][0]]
            else:
                initcond = [hf[i](0) for i in self.allvars];
            t = np.arange(t_i, t_f, dt)
            ode = ode23(self, odeint(self.eqns, initcond, t, args=(ts, qs)), t);
            return ode
        else:
            dde = dde23(eqns=self.eqnsdelay(case), params=self.paramsnd(ts, qs))
            dde.set_sim_params(tfinal=t_f, dtmax=dt, AbsTol=10 ** -6, RelTol=10 ** -3)
            if (len(hf) > 0):
                dde.hist_from_funcs(hf, 50)
            else:
                dde.hist_from_funcs(self.histfunc(case, qs), 50)
            dde.run()
            return dde

    def Jacobian(self, var):
        m, d, DNA0, DNA1, DNA2 = var;
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
        return J;

    #TODO return eigvals of all steady states
    def eigenvalues(self, ax=0, color=['k', 'b', 'g', 'r'], number=0):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        if (max(self.ss['m']) > 1e30):
            return [np.nan]
        for i in range(len(self.ss['m'])):
            J = self.Jacobian([self.ss['m'][i], self.ss['d'][i], self.ss['DNA0'][i], self.ss['DNA1'][i], self.ss['DNA2'][i]]);

            eigvals = np.linalg.eigvals(J)
            if (ax != 0):
                for eigval in eigvals[0-number:]:
                    ax.plot(eigval.real, eigval.imag, c=color[i], marker='o')
                ax.axvline(0, color='k', linestyle=':')
                ax.set_xlabel("Real")
                ax.set_ylabel("Imaginary")
        return eigvals[0-number:]

    def imageig(self):
        eig = (self.eigenvalues()).flatten()
        return np.any(eig.imag != 0)