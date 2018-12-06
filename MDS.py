from singlegenesystem import *

class MDS(SG):
    def __init__(self, params={}):
        if(len(params)==0):
            self.beta = 0
            self.gammam = 0
            self.gammamRNA = 0
            self.gammad = 0
            self.alphaass = 0
            self.alphadiss = 0
            self.Kdim = 0
            self.taum = 0
            self.taumRNA = 0
            self.DNAtot = 0
            self.phi0 = 0
            self.fm = 0
            self.fd = 0
            self.Km = 0
            self.Kd = 0
            self.kbm = 0
            self.kum = 0
            self.kbd = 0
            self.kud = 0
            self.vol = 0
        else:
            self.beta = params['beta'];  # translation rate
            self.gammam = params['gammam'];  # degradation rate monomer
            self.gammamRNA = params['gammamRNA'];  # degradation rate mRNA
            self.gammad = params['gammad'];  # degradation rate dimer
            self.alphaass = params['alphaass'];
            if ('alphadiss' in params):
                self.alphadiss = params['alphadiss'];
                self.Kdim = self.alphadiss / self.alphaass;
            elif ('Kdim' in params):
                self.Kdim = params['Kdim'];
                self.alphadiss = self.alphaass * self.Kdim;self.taum = params['taum'];  # time needed to translate monomer
            self.taumRNA = params['taumRNA'];  # time needed to transcribe mRNA
            self.DNAtot = params['DNAtot'];  # total amount of DNA copies
            self.phi0 = params['phi0'];  # basal transcription rate
            self.fm = params['fm'];  # activation fold when monomer is bound
            self.fd = params['fd'];  # repression fold when dimer is bound
            if('Km' in params):
                self.Km = params['Km'];  # binding constant monomer on DNA (k_on/k_off)
            else:
                self.Km = params['kbm']/params['kum']
            if('Kd' in params):
                self.Kd = params['Kd'];  # binding constant dimer on DNA (k_on/k_off)
            else:
                self.Kd = params['kbd']/params['kud']
            self.kbm = params['kbm'];  # k_on of monomer
            if('kum' in params):
                self.kum = params['kum']
            else:
                self.kum = self.kbm / self.Km;  # k_off of monomer
            self.kbd = params['kbd'];  # k_on of dimer
            if('kud' in params):
                self.kud = params['kud'];
            else:
                self.kud = self.kbd / self.Kd;  # k_off of dimer
            self.vol = params['vol'];  # volume
        self.DNAstates = np.array(['DNA0', 'DNAm', 'DNAd']);
        self.allvars = np.array(['DNA0', 'DNAm', 'mRNA', 'm', 'd']);
        self.taum = 0
        self.taumRNA = 0
        self.ss = {};  # steady state
        self.nameparameters= np.array(['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'fm','fd', 'Km', 'Kd', 'kbm', 'kbd']);

    def outRange(self):
        line = "";
        if self.beta < 10 ** -2.8:
            line += "beta low ";
        elif self.beta > 10 ** 2.8:
            line += "beta high ";
        if self.gammam < 10 ** -2.8:
            line += "gammam low ";
        elif self.gammam > 10 ** 2.8:
            line += "gammam high ";
        if self.gammamRNA < 10 ** -2.8:
            line += "gammamRNA low ";
        elif self.gammamRNA > 10 ** 4.8:
            line += "gammamRNA high ";
        if self.gammad < 10 ** -6.8:
            line += "gammad low ";
        elif self.gammad > 10 ** 2.8:
            line += "gammad high ";
        if self.alphaass < 10 ** -2.8:
            line += "alphaass low ";
        elif self.alphaass > 10 ** 2.8:
            line += "alphaass high ";
        if self.alphadiss < 10 ** -6.8:
            line += "alphadiss low ";
        elif self.alphadiss > 10 ** 2.8:
            line += "alphadiss high ";
        if self.phi0 < 10 ** -2.8:
            line += "phi0 low ";
        elif self.phi0 > 10 ** 2.8:
            line += "phi0 high ";
        if self.fm < 1:
            line += "fa low ";
        elif self.fm > 10 ** 2:
            line += "fa high ";
        if self.fd < 10 ** -3:
            line += "fr low ";
        elif self.fd > 1:
            line += "fr high ";
        if self.kbm < 10 ** -2.8:
            line += "kbm low ";
        elif self.kbm > 10 ** 2.8:
            line += "kbm high ";
        if self.kbd < 10 ** -2.8:
            line += "kbd low ";
        elif self.kbd > 10 ** 2.8:
            line += "kbd high ";
        if self.kum < 10 ** -2.8:
            line += "kum low ";
        elif self.kum > 10 ** 5.8:
            line += "kum high ";
        if self.kud < 10 ** -2.8:
            line += "kud low ";
        elif self.kud > 10 ** 5.8:
            line += "kud high ";
        return line;

    # return string of paramters
    def parameterline(self):
        string = "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.2E\n" % (
            self.beta, self.gammam, self.gammamRNA, self.gammad, self.alphaass, self.alphadiss, self.phi0, self.fm,
            self.fd, self.Km, self.Kd, self.kbm, self.kbd);
        return string;

    # change the parameters to new values
    def setparameters(self, params):
        self.beta = params['beta']
        self.gammam = params['gammam'];
        self.gammamRNA = params['gammamRNA'];
        self.gammad = params['gammad'];
        self.alphaass = params['alphaass'];
        self.alphadiss = params['alphadiss'];
        self.taum = params['taum'];
        self.taumRNA = params['taumRNA'];
        self.DNAtot = params['DNAtot'];
        self.phi0 = params['phi0'];
        self.fm = params['fm'];
        self.fd = params['fd']
        if ('Km' in params):
            self.Km = params['Km'];  # binding constant monomer on DNA (k_on/k_off)
        else:
            self.Km = params['kbm'] / params['kum']
        if ('Kd' in params):
            self.Kd = params['Kd'];  # binding constant dimer on DNA (k_on/k_off)
        else:
            self.Kd = params['kbd'] / params['kud']
        self.kbm = params['kbm'];  # k_on of monomer
        if ('kum' in params):
            self.kum = params['kum']
        else:
            self.kum = self.kbm / self.Km;  # k_off of monomer
        self.kbd = params['kbd'];  # k_on of dimer
        if ('kud' in params):
            self.kud = params['kud'];
        else:
            self.kud = self.kbd / self.Kd;  # k_off of dimer
        self.vol = params['vol'];
        self.ss.clear();

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
        self.fm = parent.fm
        self.fd = parent.fd
        self.kbm = parent.kbm
        self.kbm = parent.kbm
        self.kud = parent.kum
        self.kud = parent.kud
        self.Km = self.kbm / self.kum;
        self.Kd = self.kbd / self.kud;

    # set parameters to mutant of parent
    def smallmutant(self, parent, constraints={}, SsLrpB=False):
        idx = int(rd.uniform(0, 14))
        mutation = np.ones(14);
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
            self.kbm = parent.kbm * mutation[9];
            self.kum = parent.kum * mutation[10];
            self.kbd = parent.kbd * mutation[11];
            self.kud = parent.kud * mutation[12];
            self.Km = self.kbm / self.kum;
            self.Kd = self.kbd / self.kud;
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
            self.fm = max(min(parent.fm * mutation[7], constraints['f'][1]), constraints['f'][0]);#todo
            self.fd = max(min(parent.fd * mutation[8], constraints['f'][1]), constraints['f'][0]);#todo
            #self.f12 = max(min(parent.f12 * mutation[9], max(self.f1, self.f2)), constraints['f'][0]);#todo
            self.kbm = max(min(parent.kbm * mutation[10], constraints['kb'][1]), constraints['kb'][0]);
            self.kum = max(min(parent.kum * mutation[11], constraints['ku'][1]), constraints['ku'][0]);
            self.kbd = max(min(parent.kbd * mutation[12], constraints['kb'][1]), constraints['kb'][0]);
            self.kud = max(min(parent.kud * mutation[13], constraints['ku'][1]), constraints['ku'][0]);
            self.Km = self.kbm / self.kum;
            self.Kd = self.kbd / self.kud;
            self.vol = parent.vol;
        self.ss.clear();
        return

    def largemutant(self, parent, constraints={}, SsLrpB=False):
        idx = int(rd.uniform(0, 14))
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
            self.fm = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo
        elif (idx == 8):
            self.fd = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1])))#todo
        elif (idx == 10):
            self.kbm = np.exp(rd.uniform(np.log(constraints['kb'][0]), np.log(constraints['kb'][1])))
        elif (idx == 11):
            self.kbd = np.exp(rd.uniform(np.log(constraints['kb'][0]), np.log(constraints['kb'][1])))
        elif (idx == 12):
            self.kum = np.exp(rd.uniform(np.log(constraints['ku'][0]), np.log(constraints['ku'][1])))
        elif (idx == 13):
            self.kud = np.exp(rd.uniform(np.log(constraints['ku'][0]), np.log(constraints['ku'][1])))
        self.Km = self.kbm/self.kum
        self.Kd = self.kbd/self.kud
        self.ss.clear();
        return

    # calculate and return the steady state for the parameters
    def setsteadystate(self):
        # expression of dimer, DNA and mRNA as function of monomer
        fd = lambda x: self.alphaass * x ** 2 / (self.alphadiss + self.gammad);
        fDNA0 = lambda x: self.DNAtot / (self.Kd * fd(x) + self.Km * x + 1);
        fDNAm = lambda x: self.Km * x * fDNA0(x);
        fmRNA = lambda x: self.phi0 * (
            (1 - self.fd) * fDNA0(x) + (self.fm - self.fd) * fDNAm(x) + self.fd * self.DNAtot) / self.gammamRNA

        # dmRNA/dt = 0 as a function of monomer

        func = lambda x: ((self.phi0 * (1 - self.fd)) * (self.DNAtot / (self.Kd * fd(x) + self.Km * x + 1)) + (
            self.phi0 * (self.fm - self.fd)) * (self.DNAtot * self.Km * x / (
            self.Kd * fd(x) + self.Km * x + 1)) + self.phi0 * self.fd * self.DNAtot - self.gammamRNA / self.beta * (
                              self.gammam * x) - self.gammamRNA / self.beta * 2 * (
                              self.alphaass * x ** 2 - self.alphadiss * fd(x)))

        if (False):
            x = np.logspace(-7, 3, 2000)
            #x = np.linspace(0, 100, 2000)
            plt.figure('%.3E %.3E %.3E'%(self.fm, self.fd))
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
            N = 5000
            x = np.logspace(-25, 25, N);
            zeroes = np.where(np.diff(np.sign(func(x))))[0];
            if (len(zeroes) > 1):  # multistability
                zeroes = np.append(np.append([0], zeroes), [N - 1])
                m = np.array([brentq(func, (x[zeroes[i]] + x[zeroes[i - 1] + 1]) / 2,
                                     (x[zeroes[i + 1]] + x[zeroes[i] + 1]) / 2, xtol=1e-10, rtol=1e-5, maxiter=200) for
                              i in range(1, len(zeroes) - 1)])
            else:
                m = np.array([brentq(func, 1e-25, 1e25, xtol=1e-10, rtol=1e-5, maxiter=200)])

        d = fd(m);
        DNA0 = fDNA0(m);
        DNAm = fDNAm(m);
        mRNA = fmRNA(m);
        self.ss = {'DNA0': DNA0, 'DNAm': DNAm, 'DNAd': self.DNAtot - DNA0 - DNAm, 'mRNA': mRNA, 'm': m, 'd': d}
        return self.ss

    # differential equations
    def eqns(self, var, t, ts=1.0, qs=1.0, fastdimer=False):
        correct_for_zero = False;
        if (not fastdimer):
            DNA0, DNAm, mRNA, m, d = var;
            if (correct_for_zero):
                if (DNA0 < 0): DNA0 = 0; print("DNA0<0");
                if (DNAm < 0): DNAm = 0; print("DNAm<0");
                if (DNAm + DNA0 > self.DNAtot + 1e-10):
                    print("DNAm+DNA0>DNAtot");
                    DNAm *= self.DNAtot / (DNAm + DNA0);
                    DNA0 = self.DNAtot - DNAm;
                if (mRNA < 0): mRNA = 0; print("mRNA<0");
                if (m < 0): m = 0; print("m<0");
                if (d < 0): d = 0; print("d<0");

            dDNA0 = -(self.kbm * qs * ts * m + self.kbd * qs * ts * d) * DNA0 + self.kum * ts * DNAm + self.kud * ts * (
                self.DNAtot / qs - DNA0 - DNAm)
            dDNAm = self.kbm * ts * qs * m * DNA0 - self.kum * ts * DNAm;
            dmRNA = self.phi0 * ts * (
                DNA0 + self.fm * DNAm + self.fd * (self.DNAtot / qs - DNA0 - DNAm)) - self.gammamRNA * ts * mRNA
            dm = -self.kbm * ts * qs * m * DNA0 + self.kum * ts * DNAm + self.beta * ts * mRNA - 2 * self.alphaass * ts * qs * pow(
                m, 2) + 2 * self.alphadiss * ts * d - self.gammam * ts * m
            dd = -self.kbd * ts * qs * d * DNA0 + self.kud * ts * (
                self.DNAtot / qs - DNA0 - DNAm) + self.alphaass * ts * qs * pow(m,
                                                                                2) - self.alphadiss * ts * d - self.gammad * ts * d
            return [dDNA0, dDNAm, dmRNA, dm, dd];
        else:
            mRNA, m = var;
            d = self.alphaass * m ** 2 / (self.alphadiss + self.gammad);
            DNA0 = (-self.Kd * d + self.DNAtot) / (1 + self.Km * m);
            DNAm = (-self.Kd * d + self.DNAtot) / (1 + 1 / (self.Km * m));
            dmRNA = self.phi0 * DNA0 + self.phi0 * self.fm * DNAm + self.phi0 * self.fd * (
                self.DNAtot - DNA0 - DNAm) - self.gammamRNA * mRNA
            dm = - self.kbm * m * DNA0 + self.kum * DNAm + self.beta * mRNA - 2 * self.alphaass * pow(m,
                                                                                                      2) + 2 * self.alphadiss * d - self.gammam * m
            return [dmRNA, dm];

    def reactionmatrix(self):
        reactions = np.matrix([[-1, -1, 1, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0], [-1, 0, 1, 0, 0, 0, 1, -2, 2, -1, 0],
                               [0, -1, 0, 1, 0, 0, 0, 1, -1, 0, -1]])
        return reactions

    def propensities(self, var):
        DNA0, DNAm, mRNA, m, d = var;
        DNAd = self.DNAtot - DNA0 - DNAm

        prop = np.array([self.kbm * DNA0 * m, self.kbd * d * DNA0, self.kum * DNAm, self.kud * DNAd,
                         self.phi0 * (DNA0 + self.fm * DNAm + self.fd * DNAd), self.gammamRNA * mRNA,
                         self.beta * mRNA, self.alphaass * m * (m - 1), self.alphadiss * d, self.gammam * m,
                         self.gammad * d])
        if (np.any(prop < 0)):
            print(prop)
        return prop

    # delay differential equations in dictionary format for pydelay
    def eqnsdelay(self, fastdimer=False): #TODO
        if (not fastdimer):
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

    # history/initial conditions of variables
    def histfunc(self, qs, fastdimer=False):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        if (not fastdimer):
            hf = {'m': lambda t: self.ss['m'] / qs * 1.2, 'd': lambda t: self.ss['d'] / qs * 1.2,
                  'mRNA': lambda t: self.ss['mRNA'] / qs * 1.2, 'DNA0': lambda t: self.ss['DNA0'] / qs,
                  'DNAm': lambda t: self.ss['DNAm'] / qs, }
        else:
            hf = {'m': lambda t: self.ss['m'] / qs * 1.1, 'mRNA': lambda t: self.ss['mRNA'] / qs * 0.9}
        return hf

    # calculate quasi steady state with input monomer or dimer
    def quasisteadystate(self, x, monomer=True):
        if (monomer):
            m = x;
            d = self.alphaass * m ** 2 / (self.alphadiss + self.gammad);
        else:
            d = x;
            m = np.sqrt((self.alphadiss + self.gammad) * d / self.alphaass);

        # calculate DNA concentrations in quasi steady state
        DNA0 = self.DNAtot / (self.Kd * d + self.Km * m + 1);
        DNAm = self.Km * m * DNA0;
        DNAd = self.Kd * d * DNA0;

        qss = {'DNA0': DNA0, 'DNAm': DNAm, 'DNAd': DNAd, 'd': d, 'm': m};
        return qss

    # transcription rate in quasi steady state with input monomer or dimer
    def transcriptionrate(self, x, monomer=True):
        ssDNA = self.quasisteadystate(x, monomer);
        if (not len(self.ss) > 0):
            self.setsteadystate();
        return self.phi0 * ssDNA['DNA0'] + self.phi0 * self.fm * ssDNA['DNAm'] + self.phi0 * self.fd * ssDNA['DNAd']

    def translationrate(self, x, monomer=True):
        return self.beta * self.transcriptionrate(x, monomer) / self.gammamRNA

    def paramsnd(self, ts, qs):
        par = {  # nondimensionalized
            'beta': self.beta * ts, 'gammam': self.gammam * ts, 'gammamRNA': self.gammamRNA * ts,
            'gammad': self.gammad * ts, 'alphaass': self.alphaass * ts * qs, 'alphadiss': self.alphadiss * ts,
            'taum': self.taum / ts, 'taumRNA': self.taumRNA / ts, 'DNAtot': self.DNAtot / qs, 'phi0': self.phi0 * ts,
            'fm': self.fm, 'fd': self.fd, 'Km': self.Km * qs, 'Kd': self.Kd * qs, 'kbm': self.kbm * qs * ts,
            'kum': self.kum * ts, 'kbd': self.kbd * qs * ts, 'kud': self.kud * ts, 'volume': self.vol, 'ts': ts,
            'qs': qs}
        return par

    def nullClineM(self, x, monomer=True):  # TODO: IF NOT MONOMER RETURN DIMER STATE???
        qss = self.quasisteadystate(x, monomer);
        # values of mRNA for which dm/dt are 0 :
        m0 = (self.kbm * qss['m'] * qss['DNA0'] - self.kum * qss['DNAm'] + 2 * self.alphaass * qss[
            'm'] ** 2 - 2 * self.alphadiss * qss['d'] + self.gammam * qss['m']) / self.beta;
        return m0;

    def nullClineMRNA(self, x, monomer=True):
        qss = self.quasisteadystate(x, monomer);
        # values of mRNA for which dmRNA/dt are 0 :
        mRNA0 = self.phi0 * (qss['DNA0'] + self.fm * qss['DNAm'] + self.fd * qss['DNAd']) / self.gammamRNA;
        return mRNA0;

    def timelapse(self, t_f, dt, t_i=0, hf={}, fastdimer=False):
        ts = 1.0;
        qs = 1.0;  # SCALING
        if (self.taum == 0 and self.taumRNA == 0):
            if (len(hf) == 0):  # no initial conditions given, start from perturbation of steady state
                if (not len(self.ss) > 0):
                    self.setsteadystate();
                if (not fastdimer):
                    y0 = [1, 0, 0, 1, 0];
                    #[self.ss['DNA0'], self.ss['DNAm'], self.ss['mRNA'] * 1.05, self.ss['m'] * 1.05,
                         # self.ss['d'] * 1.05]
                else:
                    y0 = [self.ss['mRNA'] * 1.2, self.ss['m'] * 1.2]
            else:
                if (not fastdimer):
                    y0 = [hf['DNA0'](0), hf['DNAm'](0), hf['mRNA'](0), hf['m'](0), hf['d'](0)]
                else:
                    y0 = [hf['mRNA'](0), hf['m'](0)]
            t = np.arange(t_i, t_f, dt);
            ode = ode23(self, odeint(self.eqns, y0, t, args=(ts, qs, fastdimer)), t);
            return ode;
        else:
            dde = dde23(eqns=self.eqnsdelay(), params=self.paramsnd(ts, qs))
            dde.set_sim_params(tfinal=t_f, dtmax=dt, AbsTol=10 ** -6, RelTol=10 ** -3)
            if (len(hf) > 0):
                dde.hist_from_funcs(hf, 50)
            else:
                dde.hist_from_funcs(self.histfunc(qs), 50)
            dde.run()
            return dde

    def Jacobian(self, var, fastdimer=False):
        m, d, DNA0 = var;
        if (not fastdimer):
            J = [
                [-(self.kbm * m + self.kbd * d + self.kud), self.kum - self.kud, 0, -self.kbm * DNA0, -self.kbd * DNA0],
                [self.kbm * m, -self.kum, 0, self.kbm * DNA0, 0],
                [self.phi0 * (1 - self.fd), self.phi0 * (self.fm - self.fd), -self.gammamRNA, 0, 0],
                [-self.kbm * m, self.kum, self.beta, -self.kbm * DNA0 - 4 * self.alphaass * m - self.gammam,
                 2 * self.alphadiss],
                [-self.kbd * d - self.kud, -self.kud, 0, 2 * self.alphaass * m,
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

    def eigenvalues(self, ax=0, color=['k', 'b', 'g', 'r'], number=0):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        if (max(self.ss['m']) > 1e30):
            return [np.nan]
        for i in range(len(self.ss['m'])):
            J = self.Jacobian([self.ss['m'][i], self.ss['d'][i], self.ss['DNA0'][i]]);

            eigvals = np.linalg.eigvals(J)
            if (ax != 0):
                for eigval in eigvals[0-number:]:
                    ax.plot(eigval.real, eigval.imag, c=color[i], marker='o')
                ax.axvline(0, color='k', linestyle=':')
                ax.set_xlabel("Real")
                ax.set_ylabel("Imaginary")
        return eigvals[0-number:]