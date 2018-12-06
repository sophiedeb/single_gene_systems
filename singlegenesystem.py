import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from scipy.optimize import brentq # for daughter classes
from scipy.integrate import odeint # for daughter classes

# TODO qs ts delay oscillationParams in DDO and MDO as in DDDO

class ode23:
    def __init__(self, oscillator, conc, t, reaccount=[]):
        conc[conc < 0.0] = 0.0;  # concentrations cannot be negative, integration errors can become negative
        self.t = t;
        self.reaccount = reaccount;
        if (oscillator.__class__.__name__ ==  'DDS'):
            if (len(conc[0]) == 2):  # fastdimer
                self.sol = {'mRNA': conc[:, 0], 'm': conc[:, 1]}
            else:
                self.sol = {'DNA0': conc[:, 0], 'DNA1': conc[:, 1], 'DNA2': conc[:, 2],
                            'DNA12' : oscillator.DNAtot - conc[:,0] - conc[:,1] - conc[:,2],
                            'mRNA': conc[:, 3], 'm': conc[:, 4], 'd': conc[:, 5]}
        elif (oscillator.__class__.__name__ == 'DS'):
            self.sol = {'DNA0': conc[:, 0], 'mRNA': conc[:, 1], 'm': conc[:, 2], 'd': conc[:, 3]}
        elif (oscillator.__class__.__name__ == 'MDS'):
            if (len(conc[0]) == 2):  # fastdimer
                self.sol = {'mRNA': conc[:, 0], 'm': conc[:, 1]}
            else:
                self.sol = {'DNA0': conc[:, 0], 'DNAm': conc[:, 1], 'mRNA': conc[:, 2], 'm': conc[:, 3],
                            'd': conc[:, 4], 'DNAd' : oscillator.DNAtot - conc[:,0] - conc[:,1]}
        elif (oscillator.__class__.__name__ == 'DDDS'):
            self.sol = {'DNA0': conc[:, 0], 'DNA1': conc[:, 1], 'DNA2': conc[:, 2], 'DNA3': conc[:, 3],
                        'DNA12': conc[:, 4], 'DNA13': conc[:, 5], 'DNA23': conc[:, 6],
                        'DNA123' : oscillator.DNAtot - conc[:,0] - conc[:,1] - conc[:,2] - conc[:,3] - conc[:,4] - conc[:,5] - conc[:,6],
                        'mRNA': conc[:, 7], 'm': conc[:, 8], 'd': conc[:, 9]}
        elif isinstance(oscillator, MS):
            self.sol = {'DNA0': conc[:, 0], 'mRNA': conc[:, 1], 'm': conc[:, 2]}

class SG:
    def __init__(self):
        self.beta = None;
        self.gammam = None;
        self.gammamRNA = None;
        self.gammad = None;
        self.alphaass = None;
        self.alphadiss = None;
        self.taum = None;
        self.taumRNA = None;
        self.DNAtot = None;
        self.vol = params['vol'];
        self.ss = {};
        self.DNAtot = 0;
        self.DNAstates = [];
        self.allvars = [];
        return

    def nullClineM(self, x, monomer=True):
        return

    def nullClineMRNA(self, x, monomer=True):
        return

    def mutant(self, parent):
        return

    def setsteadystate(self):
        return

    def setparameters(self, params):
        return;

    def steadystate(self):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        return self.ss

    def eqns(self, var, t=0.0, ts=1.0, qs=1.0):
        return;

    def histfunc(self, qs, fastdimer):
        return;

    def quasisteadystate(self, x, monomer=True):
        return;

    def transcriptionrate(self, m):
        return;

    def paramsnd(self, ts, qs):
        return;

    def extrema(self, t_f, dt, vars, hf={}):
        if (len(hf) == 0):
            hf = self.histfunc();
        ode = self.timelapse(t_f, dt, hf=hf).sol
        # function to remove adjacent duplicates
        remdup = lambda x: np.array([n for i, n in enumerate(x) if i == 0 or n != x[i - 1]]);
        ext = {};
        for var in vars:
            if (var in self.allvars):
                # idxvar = self.allvars.index(var);
                odevar = ode[var];
            elif (var in self.DNAstates):
                odevar = self.DNAtot - sum(
                    [ode[i] for i in self.DNAstates[:-1]]);  # ode[:-1, 0] - ode[:-1, 1] - ode[:-1, 2];
            else:
                print('variable %s not found'%var)
                return
            odevar = remdup(odevar);
            ode0 = odevar[:-1];
            ode1 = odevar[1:];
            pos = [(np.r_[False, ode1 < ode0] & np.r_[ode0 < ode1, False]) | (
                np.r_[False, ode1 > ode0] & np.r_[ode0 > ode1, False])];  # first and last elements always false
            ext[var] = ode0[pos[0][:-1]];
        return ext;

    def scores(self, t_f, dt, hf={}):
        scores = [20 for _ in range(len(self.allvars)+1)];
        if(False):
            for i in range(len(self.DNAstates)):
                scores[i] = self.scorefunction(t_f, dt, self.DNAstates[i], hf);
            for i in range(len(self.allvars) - (len(self.DNAstates)-1)):
                scores[len(self.DNAstates)+i] = self.scorefunction(t_f, dt, self.allvars[len(self.DNAstates)-1+i], hf);
        else:
            vars = self.DNAstates;
            for i in range(len(self.DNAstates)-1, len(self.allvars)):
                vars = np.append(vars, self.allvars[i])
            scores = self.scorefunction(t_f, dt, vars, hf)
        return scores;

    def scoreline(self, t_f, dt, hf={}):
        scores = self.scores(t_f, dt, hf);
        line = "";
        for score in scores:
            line += "%.3E\t"%score
        line += "\n"
        return line;

    def scorefunction(self, t_f, dt, var, hf={}):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        if (type(var) == str): #TODO str version not same as list version
            ext = self.extrema(t_f, dt, [var], hf)[var];
            S = 20;
            S -= 1/(1+np.exp((self.ss['mRNA'][0]-20.0)));
            S -= np.tanh(self.ss['d'][0]/20.0);
            S -= np.tanh((self.ss['d'][0]-self.ss['m'][0])/20.0);
            S -= np.tanh((self.ss['d'][0]-self.ss['mRNA'][0])/20.0);
            if (len(ext) >= 2):
                S -= 0.5;
            if(len(ext)>=5):
                S -= 0.5;
            if(len(ext)>=7):
                S -= 0.5
            if (len(ext) >= 9):
                S -= 0.5
            if (len(ext) > 2):
                if (len(ext) <= 11):
                    begin = 1;
                else:
                    begin = 1 #len(ext) - 11; todo
                for i in range(begin, min(11, len(ext)-1)):#len(ext) - 1):
                    S = S - abs(ext[i] - ext[i + 1]) / (ext[i] + ext[i + 1]) * min(1., abs(ext[i] - ext[i + 1]));
        else:
            S = [20 for _ in range(len(var))];
            extrema = self.extrema(t_f, dt, var, hf);
            for k in range(len(var)):
                ext = extrema[var[k]];
                S[k] -= 1 / (1 + np.exp((self.ss['mRNA'][0] - 20.0)));
                #if(isinstance(self, MDS)):
                #    S[k] -= np.tanh(self.ss['d'][0] / 20.0);
                #    S[k] -= np.tanh(self.ss['m'][0] / 20.0);
                #else:
                S[k] -= np.tanh(self.ss['d'][0] / 50.0);
                S[k] -= np.tanh((self.ss['d'][0] - self.ss['m'][0]) / 20.0);
                if (len(ext) >= 3):
                    S[k] -= 0.5;
                if (len(ext) >= 5):
                    S[k] -= 0.5;
                if (len(ext) >= 7):
                    S[k] -= 0.5
                if (len(ext) >= 9):
                    S[k] -= 0.5
                if (len(ext) > 2):
                    if (len(ext) <= 11):
                        begin = 1;
                    else:
                        begin = 1 #len(ext) - 10;
                    for i in range(begin, min(11,len(ext) - 1)):
                        S[k] = S[k] - abs(ext[i] - ext[i + 1]) / (ext[i] + ext[i + 1]) * min(1., abs(
                            ext[i] - ext[i + 1]));
                        if (ext[i] + ext[i + 1] == 0):
                            S[k] = 20.0;
            #S = np.mean(S)
        return S;

    def timelapse(self, t_f, dt, t_i=0, hf={}):
        return;

    def Jacobian(self, var):
        return;

    def eigenvalues(self, ax=0, color='k'):
        return;

    def parameterline(self):
        return

    def plotSteadystate(self, ax, monomer=True, maxx=0):
        if (not len(self.ss) > 0):
            self.setsteadystate();

        if (monomer and maxx>0):
            x = np.linspace(0, maxx, 1000)
            m = x
        elif (monomer):
            x = np.linspace(0, 5 * max(self.ss['m']), 1000)
            m = x
        elif(maxx>0):
            x = np.linspace(0, maxx, 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);
        else:
            x = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);

        ssDNA = self.quasisteadystate(m);
        for state in self.DNAstates:
            ax.plot(x, ssDNA[state] / self.DNAtot, label=state);
        ax.legend();
        ax.set_ylabel('Ratios at steady-state');
        if(monomer):
            ax.set_xlabel('Monomer concentration');
            for ssm in self.ss['m']:
                ax.axvline(ssm, color='k', linestyle=':');
        else:
            ax.set_xlabel('Dimer concentration');
            for ssd in self.ss['d']:
                ax.axvline(ssd, color='k', linestyle=':');

    def plotNullclines(self, ax, monomer=True):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        if (monomer):
            x = np.linspace(0, 20 * max(self.ss['m']), 1000);
        else:
            x = np.linspace(0, 20 * max(self.ss['d']), 1000);

        x0 = self.nullClineM(x, monomer);
        mRNA0 = self.nullClineMRNA(x, monomer);

        if (monomer):
            #x /= min(self.ss['m']);
            x0label = 'dm/dt = 0';
            xlabel = r'Monomer copy number $m$' # [%.2E M]' % min(self.ss['m']);
            if(False):
                t = [True if i%100==0 else False for i in range(1000)]
                for m in np.linspace(0,10,10):#max(x),20):
                    for mRNA in np.linspace(0,6,10):#3*max([max(x0),max(mRNA0)]),10):
                        qss = self.quasisteadystate(m, monomer)
                        # TODO only for DDS
                        _, _, _, dmRNA, dm, _ = self.eqns([qss['DNA0'], qss['DNA1'], qss['DNA2'], mRNA, m, qss['d']]);
                        ax.quiver(m, mRNA, dm, dmRNA, color='lightgrey', linewidths=0.1)
        else:
            #x /= min(self.ss['d']);
            x0label = 'dm/dt = 0';  # TODO : label
            xlabel = r'Dimer copy number' #[%.2E M]' % min(self.ss['d']);
            if(False):
                for d in np.linspace(0,60, 10): #max(x),20):
                    for mRNA in np.linspace(0,7.5, 10): #3*max([max(x0),max(mRNA0)]),10):
                        qss = self.quasisteadystate(d, monomer)
                        # TODO only for DDS
                        _, _, _, dmRNA, _, dd = self.eqns([qss['DNA0'], qss['DNA1'], qss['DNA2'], mRNA, qss['m'], d]);
                        ax.quiver(d, mRNA, dd, dmRNA, lw=0.5) #, scale=0.0001)

        #ax.plot(x, x0 / min(self.ss['mRNA']), label=x0label, color='r');  # TODO : label
        #ax.plot(x, mRNA0 / min(self.ss['mRNA']), label='dmRNA/dt = 0', color='orange');
        ax.plot(x, x0, label=x0label, color='r');  # TODO : label
        ax.plot(x, mRNA0, label='dmRNA/dt = 0', color='orange');
        ax.legend(loc=2);
        ax.set_xlabel(xlabel);
        ax.set_ylabel(r'mRNA copy number') #[%.2E M]' % min(self.ss['mRNA']));

    def plotRateBalance(self, ax, monomer=True, maxx=0):
        if (not len(self.ss) > 0):
            self.setsteadystate();

        if(monomer and maxx>0):
            x = np.linspace(0, maxx, 1000)
            m = x
        elif (monomer):
            x = np.linspace(0, 5 * max(self.ss['m']), 1000)
            m = x
        elif(maxx>0):
            x = np.linspace(0, maxx, 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);
        else:
            x = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);

        ax.plot(x, self.translationrate(m, monomer=True), label='translation');
        ax.plot(x, self.gammam * m, label='death');
        if (self.alphaass != None):
            #ax.plot(x, self.translationrate(m, monomer=True) + 2 * self.alphadiss * self.alphaass * m ** 2 / (self.alphadiss + self.gammad), label='dimerization')
            #ax.plot(x, self.gammam * m + 2 * self.alphaass * m ** 2, label='dimerization')
            ax.plot(x, 2 * self.gammad * self.alphaass / (self.alphadiss + self.gammad) * m ** 2, label='dimerization')
        ax.plot(x, self.translationrate(m, monomer=True) - self.gammam*m - 2 * self.gammad * self.alphaass / (self.alphadiss + self.gammad) * m ** 2, label='sum')

        ax.legend()
        ax.set_ylabel('Rates')
        if(monomer):
            ax.set_xlabel('Monomer concentration')
            for ssm in self.ss['m']:
                ax.axvline(ssm, color='k', linestyle=':')
        else:
            ax.set_xlabel('Dimer concentration')
            for ssd in self.ss['d']:
                ax.axvline(ssd, color='k', linestyle=':')

    def plotResponseCurve(self, ax, monomer=True, maxx=0):
        if (not len(self.ss) > 0):
            self.setsteadystate();

        if(monomer and maxx>0):
            x = np.arange(0, maxx)
            #m = x
        elif (monomer):
            x = np.linspace(0, 5 * max(self.ss['m']), 1000)
            #m = x
        elif(maxx>0):
            x = np.arange(0, maxx)
            #m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);
        else:
            x = np.arange(0, 2 * max(self.ss['d']))
            #m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);

        ax.plot(x, self.translationrate(x, monomer));

    def scorebistab(self):
        if (not len(self.ss) > 0):
            self.setsteadystate();

        #score for bistability on 10:
        if(len(self.ss['m'])<3):
            x = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * x / self.alphaass);

            x = self.translationrate(m, monomer=True) - self.gammam*m - 2 * self.gammad * self.alphaass / (self.alphadiss + self.gammad) * m ** 2

            lm = np.r_[False, x[1:] < x[:-1]] & np.r_[x[:-1] < x[1:], False] #local minima
            if(sum(lm)==0):
                score = 0;
            elif(sum(lm)>1):
                print("sum bigger than 1")
                score = 0
            else:
                score = 10.0/(1+np.exp(-2*(x[lm][0]-10)))
        else:
            scorethreshold = np.tanh((self.ss['d'][1]-self.ss['d'][0])/20.0)*np.tanh((self.ss['d'][2]-self.ss['d'][1])/20.0)
            scoremax = np.tanh(self.ss['d'][2]/100.0)*1/(1+np.exp(0.01*(self.ss['d'][2]-1500.0)))

            scoremRNAmax = 1/(1+np.exp((self.ss['mRNA'][2]-20.0)))

            d = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * d / self.alphaass);

            x = self.translationrate(m, monomer=True) - self.gammam * m - 2 * self.gammad * self.alphaass / (
            self.alphadiss + self.gammad) * m ** 2;

            posx = x[np.logical_and(d > self.ss['d'][1], d < self.ss['d'][2])]

            scoretl = np.tanh(sum(posx)/len(posx)/5.0)

            #bigx = x[np.logical_and(d > self.ss['d'][2], d < 1.15*self.ss['d'][2])]
            #scoretl2 = np.tanh(-sum(bigx)/len(bigx)/5.0)

            #_,_,_,_,it = self.parambistab()

            #scoreit = 1/(1+np.exp(0.01*(it-500.0)))

            #lm = np.r_[False, x[1:] > x[:-1]] & np.r_[x[:-1] > x[1:], False]  & np.r_[d > self.ss['d'][1]] & np.r_[d < self.ss['d'][2]] # local maximum between unstable and stable steady state
            #if(len(x[lm])==1):
            #    scoretl = np.tanh(x[lm][0]/5.0);
            #else:
            #    scoretl = 1
            #    print("multiple maxima")
            score = 10 + 2.5*scorethreshold + 2.5*scoremax + 2.5*scoremRNAmax + 2.5*scoretl #+ 2.5*scoretl + 2.5*scoretl2;
        return 20  - score;

    def parambistab(self, stoch = False): #TODO changed definition of induction time, changed return value
        if (not len(self.ss) > 0):
            self.setsteadystate();

        if(len(self.ss['d'])>2):
            d = np.linspace(0, 2 * max(self.ss['d']), 1000)
            m = np.sqrt((self.alphadiss + self.gammad) * d / self.alphaass);

            x = self.translationrate(m, monomer=True) - self.gammam * m - 2 * self.gammad * self.alphaass / (
                self.alphadiss + self.gammad) * m ** 2;

            posx = x[np.logical_and(d > self.ss['d'][1], d < self.ss['d'][2])]
            bigx = x[np.logical_and(d > self.ss['d'][2], d < 1.15 * self.ss['d'][2])]

            t_f = 10000;
            dt = 1E-2;

            hf = {}

            def set_lambda(par):  # need this function otherwise problems with loops over dictionaries of lambda functions
                if par=="d":
                    return lambda t: 1.1 * self.steadystate()[par][1]
                else:
                    return lambda t: self.steadystate()[par][1]
            for par in self.allvars:
                hf[par] = set_lambda(par)

            ts = self.timelapse(t_f, dt, hf=hf)
            tsd = ts.sol['d']

            #fig = plt.figure()
            #ax = fig.add_subplot(111)
            #self.plottimelapse(ax, ts, t_f, dt)
            #plt.show()

            it = len(tsd[tsd<0.90*self.steadystate()['d'][2]])*dt
            if(False):
                for par in self.allvars:
                    print(par, hf[par](0))
                print(tsd)
                fig = plt.figure("ts %.3E"%it)
                ax = fig.add_subplot(1,1,1)
                self.plottimelapse(ax, ts, t_f, dt)
                ax.axhline(0.95*self.steadystate()['d'][2])
                plt.show()

            if(stoch):
                itstoch = np.zeros(5)
                for i in range(5):
                    itstoch[i] = self.stochasticinductiontime()

                return it, itstoch[0], itstoch[1], itstoch[2], itstoch[3], itstoch[4]
            else:
                return it
        else:
            return np.nan

    def plottimelapse(self, ax, de, t_f, dt, t_i=0, DNA=False, axncm=0, axncd=0, legend=True):
        if (self.taum == 0 and self.taumRNA == 0):
            t = de.t;
            sol = de.sol;
        else:  # use pydelay
            if (dt < (t_f - t_i) / 5000):  # resample for plot if interval is too small
                dt = (t_f - t_i) / 5000;
            sol = de.sample(t_i, t_f, dt);
            t = np.arange(t_i, t_f, dt);

        qs = 1.0;  # de.params['qs'];
        ts = 1.0;  # de.params['ts'];

        if (not DNA):
            m = sol['m']
            mRNA = sol['mRNA']
            if (self.alphaass != None):
                if ('d' in sol):
                    d = sol['d'];
                else:
                    d = self.quasisteadystate(m)['d'];
                conc = [m, d, mRNA];
                labels = ['Monomer', 'Dimer', 'mRNA'];
                sslabel = ['m', 'd', 'mRNA'];
            else:
                conc = [m, mRNA];
                labels = ['Monomer', 'mRNA'];
                sslabel = ['m', 'mRNA'];
        elif ('DNA0' in sol):
            ax.set_ylim([0, self.DNAtot*1.1]);
            conc = [sol[state] for state in self.DNAstates[:-1]];
            conc += [self.DNAtot - sum(conc)];
            labels = self.DNAstates;
            sslabel = labels;
        else:
            ax.set_ylim([0, self.DNAtot*1.1]);
            m = sol['m']
            mRNA = sol['mRNA']
            qss = self.quasisteadystate(m);
            conc = [qss[state] for state in self.DNAstates];
            labels = self.DNAstates;
            sslabel = labels;
        if (ax != 0):
            for i in [0]: #range(len(conc)): #[0]
                l = ax.plot(t * ts, conc[i], label=labels[i], c='k');
                if (not len(self.ss) > 0):
                    self.setsteadystate();
                    print('change in plot')
                if (sslabel[i] in self.ss):
                    ls = ['--', ':','-.']
                    for j in range(len(self.ss[sslabel[i]])):
                        ax.axhline(self.ss[sslabel[i]][j] * qs, color=l[0].get_color(), linestyle=ls[j], linewidth=0.8)
            ax.set_xlabel('Time [min]') #, fontsize=12);
            if(not DNA):
                ax.set_ylabel('Copy number') #, fontsize=12, rotation=0);
                #ax.yaxis.set_label_coords(-0.02, 1.05)
            # ax.set_ylabel('Concentration [%.2EM]'%qs);
            if(legend):
                ax.legend(loc=1, fontsize=11, handlelength=1);
            # ax2 = ax.twinx();
            # mn, mx = ax.get_ylim();
            # ax2.set_ylim(mn*self.vol*qs*navo, mx*self.vol*qs*navo);
            # ax2.set_ylabel('Copy number');
            ax.tick_params(axis='both', which='major', labelsize=11)  # , **calfont)
        if (axncm != 0):
            #axncm.plot(m[len(m) / 3:] / min(self.ss['m']), mRNA[len(mRNA) / 3:] / min(self.ss['mRNA']));
            axncm.plot(m[len(m) / 3:], mRNA[len(mRNA) / 3:], label='timeseries');
            axncm.set_xlim([0, 1.2*max(m[len(m)/3:])]);
        if (axncd != 0):
            #axncd.plot(d[len(d) / 3:] / min(self.ss['d']), mRNA[len(mRNA) / 3:] / min(self.ss['mRNA']));
            axncd.plot(d[len(d) / 3:], mRNA[len(mRNA) / 3:]);
            axncd.set_xlim([0, 1.2*max(d[len(d)/3:])]);
        return ax

    def stochasticinductiontime(self): # TODO implement initial conditions
        if (not len(self.ss) > 0):
            self.setsteadystate();
        if (self.taum == 0 and self.taumRNA == 0):
            x = np.array([[1, 0, 0, 0, 0, 0, 0, int(self.ss['mRNA'][1])+1, int(self.ss['m'][1])+1,int(self.ss['d'][1])+15]])

            t = 0
            reactions = self.reactionmatrix()
            maxd = 0.95*max(self.steadystate()['d'])
            tt = 50;
            ts = np.zeros(1)
            dt = 0.1
            tt2 = 0.1
            xs = np.copy(x)

            print('start stochastic')
            while (t < 5000 and x[0][-1] < maxd):
                prop = self.propensities(x[0]);
                proptot = sum(prop);

                reac = np.random.choice(len(prop), p=prop / proptot);
                Km = np.array([i == reac for i in range(len(prop))]);
                ru = np.random.uniform();
                tau = 1 / proptot * np.log(1 / ru);

                t += tau;

                while (t >= tt2):  # save before adapting x, correct for K=1, approximation for K>1
                    ts = np.append(ts, t)
                    xs = np.vstack((xs, x))
                    tt2 += dt;
                if(t>tt):
                    print(t, x)
                    tt = 50*(int(t/50)+1)

                x += reactions.dot(Km);
            #plt.figure()
            #plt.plot(ts,xs)
            #plt.show()

            print('end stochastic', t)
            return t

    def stochastictimelapse(self, t_f, dt, t_i=0, hf={}): # TODO implement initial conditions
        if (not len(self.ss) > 0):
            self.setsteadystate();
        if (self.taum == 0 and self.taumRNA == 0):
            x = np.zeros((1,len(self.allvars)))
             # TODO implement initial conditions!
            if(len(hf)!=0):
                print('code yet to be written for initial conditions not equal to steady state')
            else:
                #x = hf.copy()
                #DNA0, DNA1, DNA2, DNA3, DNA12, DNA23, DNA13, mRNA, m, d = var;
                hf=self.ss
                for idx, par in enumerate(self.allvars):
                    if par == 'DNA0':
                        x[0,idx]=1
                    elif(par in self.DNAstates and par != 'DNA0'):
                        x[0,idx]=0
                    else:
                        x[0,idx]=int(hf[par][0])

            t = t_i
            tt = t_i + dt;
            xts = x.copy()
            reactions = self.reactionmatrix()
            reaccount = np.zeros((reactions.shape)[1])
            atH = False;
            print('bug',self.allvars)
            idxd = np.where(self.allvars=='d')[0][0]
            maxd = 0.95*max(self.steadystate()['d'])

            while (t < t_f):
                prop = self.propensities(x[0]);
                proptot = sum(prop);

                reac = np.random.choice(len(prop), p=prop / proptot);
                Km = np.array([i == reac for i in range(len(prop))]);
                ru = np.random.uniform();
                tau = 1 / proptot * np.log(1 / ru);
                if(not atH):
                    reaccount[reac] += 1

                t += tau;

                while (t >= tt and tt < t_f-dt/2):  # save before adapting x, correct for K=1, approximation for K>1
                    xts = np.vstack((xts, x))
                    tt += dt;
                    if(tt%0.01==0):
                        print(t,tau, x)
                x += reactions.dot(Km);
                if(not atH and xts[-1][idxd] > maxd):
                    atH = True;
            ode = ode23(self, xts, np.arange(t_i, t_f, dt), reaccount)
            return ode

    def plotDNATimeAveraged(self, ax, ode, timestep):
        def helper(X, t, ti, tf):
            if (ti < t[0]):
                idxi = 0;
            else:
                idxi = np.where((t[:-1] <= ti) & (t[1:] > ti))[0][0]
            if (tf >= t[-1]):
                idxf = len(t) - 1;
            elif (tf < t[0]):
                return 0
            else:
                idxf = np.where((t[:-1] <= tf) & (t[1:] > tf))[0][0] + 1
            th = t[idxi:idxf];
            th[0] = ti;
            th[-1] = tf;
            Xh = X[idxi:idxf];
            res = np.sum([(te - tb) * nX for tb, te, nX in zip(th[:-1], th[1:], Xh[:-1])])
            return res

        t = ode.t

        newt = np.arange(0, t[-1], timestep)

        DNA = {}
        DNA[self.DNAstates[-1]] = np.ones(len(t))*self.DNAtot;
        for i in self.DNAstates[:-1]:
            DNA[i] = ode.sol[i];
            DNA[self.DNAstates[-1]] -= DNA[i];
        print(DNA['DNA12'])

        TADNA = {};

        for i in self.DNAstates:
            TADNA[i] = np.zeros(len(newt))

        for i in range(1, len(newt)):
            tDNA = [helper(DNA[j], t, newt[i - 1], newt[i]) for j in self.DNAstates];
            TADNA[self.DNAstates[np.argmax(tDNA)]][i] = 1;

        for i in self.DNAstates:
            ax.fill_between(newt, 0, TADNA[i], edgecolor="none", label=i, step="post")
        ax.legend()
        ax.set_yticks([])
        ax.tick_params(axis='both', which='major', labelsize=18)  # , **calfont)
        ax.set_xlabel('Time (min)', fontsize=12)
        ax.legend(loc='upper center', markerscale=0.35, columnspacing=1, bbox_to_anchor=(0.5, -1.5), ncol=4, fontsize=12)

    def monotonic(self):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        x = np.linspace(min(self.steadystate()['m']), 5*max(self.steadystate()['m']), 2000);
        tx = self.translationrate(x);
        if((max(tx)-tx[-1])/max(tx) > 0.05):
            return False
        else:
            return True

    def transcriptionrate(self, x, monomer=False):
        return

    def monotonicitytype(self):
        xx = np.linspace(0, 1000, 1000)

        y = self.transcriptionrate(xx, monomer=False)
        dy = y[1:] - y[:-1]

        zero_crossings = np.where(np.diff(np.sign(dy)))[0]

        if (len(zero_crossings) == 1):
            if dy[zero_crossings[0]] > 0:
                nonmono = 1;
            elif dy[zero_crossings[0]] < 0:
                nonmono = 2;
        elif (len(zero_crossings) > 1):
            nonmono = 3;
        else:
            nonmono = 0;

        return nonmono

    def SsLrpBcompatible(self):
        xx = np.linspace(0, 1000, 1000)

        y = self.transcriptionrate(xx, monomer=False)

        return (self.monotonicitytype() == 1 and max(y)/y[0] > 2 and min(y)/y[0] < 0.5)

    def nonmonotonicity(self):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        x = np.linspace(min(self.ss['m']), max(self.ss['m']),100)
        maxtrans = max(self.transcriptionrate(x))
        endtrans = self.transcriptionrate(x[-1])
        #endtrans = self.transcriptionrate(max(self.ss['m']))
        #if(maxtrans-self.transcriptionrate(x[-1])>0):
        #    return 1 + (maxtrans-endtrans)/maxtrans
        #else:
        return (maxtrans-endtrans)/maxtrans

    def oscillationParameters(self, vars):  # eigenvalues have to have positive real part, otherwise nan
        if (not len(self.ss) > 0):
            self.setsteadystate();
        if (len(self.ss['m'])>1):
            return [np.nan, np.nan, np.nan]
        t_f = 3000;
        dt = 0.1;
        if (self.ss['mRNA'][0] > self.ss['m'][0] and self.ss['mRNA'][0] > self.ss['d'][0]):
            hvar = 'mRNA';
        elif (self.ss['m'][0] > self.ss['mRNA'][0] and self.ss['m'][0] > self.ss['d'][0]):
            hvar = 'm';
        elif (self.ss['d'][0] > self.ss['mRNA'][0] and self.ss['d'][0] > self.ss['m'][0]):
            hvar = 'd';
        if (not np.any(self.eigenvalues().real > 0)):
            return [np.nan] * 3
        while (True and t_f/dt < 1e6):
            de = self.timelapse(t_f, dt)
            if (False):
                fig, ax = plt.subplots()
                self.plottimelapse(ax, de, t_f, dt)
                plt.show()
            if (self.taum == 0 and self.taumRNA == 0):
                sol = de.sol;
            else:
                sol = de.sample(0, t_f, dt);
            qs = 1;
            ts = 1;
            #qs = de.params['qs'];
            #ts = de.params['ts'];
            locmax = argrelextrema(sol[hvar], np.greater)[0]
            locmin = argrelextrema(sol[hvar], np.less)[0]
            if (len(locmax) > 5 and len(locmin) > 5 and (
                        abs(sol[hvar][locmax[-4]] - sol[hvar][locmax[-3]]) < 0.01 * sol[hvar][locmax[-3]]) and (
                        abs(sol[hvar][locmin[-4]] - sol[hvar][locmin[-3]]) < 0.01 * sol[hvar][locmin[-3]]) and (
                        (2 * locmax[-3] - locmax[-2] - locmax[-4]) < 3)):  # [-1] can be due to last element in array
                if (False):
                    fig, ax = plt.subplots()
                    self.plottimelapse(ax, de, t_f, dt)
                    plt.show()
                if isinstance(vars, basestring): #Python 3: str, python 2 : basestring
                    var = vars
                    if (var != hvar):
                        locmax = argrelextrema(sol[var], np.greater)[0]
                        locmin = argrelextrema(sol[var], np.less)[0]
                    period = dt * ts * (locmax[-2] - locmax[-3]);
                    ampmin = qs * sol[var][locmin[-2]];
                    ampmax = qs * sol[var][locmax[-2]];
                    return [period, ampmin, ampmax];
                else:
                    ospars = []
                    for i,var in enumerate(vars):
                        locmax = argrelextrema(sol[var], np.greater)[0]
                        locmin = argrelextrema(sol[var], np.less)[0]
                        if(len(locmax)<3 or len(locmin)<3):
                            return [np.nan, np.nan, np.nan]
                        ampmin = qs * sol[var][locmin[-2]];
                        ampmax = qs * sol[var][locmax[-2]];
                        if (i == 0):
                            period = dt * ts * (locmax[-2] - locmax[-3]);
                            ospars += [period, ampmin, ampmax];
                        else:
                            ospars += [ampmin, ampmax];
                    return ospars
            else:
                #fig, ax = plt.subplots()
                #self.plottimelapse(ax, de, t_f, dt)
                #plt.show()
                t_f = 2 * t_f;
                dt = dt/2;
        return [np.nan, np.nan, np.nan];

    def isoscillatory(self):
        eig = (self.eigenvalues()).flatten()
        return np.any(np.logical_and(eig.imag != 0, eig.real>0))
