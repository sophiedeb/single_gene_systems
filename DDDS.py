from singlegenesystem import *

def physiologicalRange():
    #PHI0 CHANGE FOR BISTABILITY, BINDING, UNBINDING
    params = {'beta': (0.01, 120), 'gammam': (1e-3, 1), 'gammamRNA': (1e-3, 10.0), 'gammad': (5e-4, 0.1),
              'alphaass': (1e-2, 1), 'alphadiss': (1e-3, 1e3),
              'phi0': (1e-2, 1e1), 'fa': (1, 100), 'fr': (1e-3, 1), 'f': (1e-3, 100),
              'kb': (0.001, 100), 'ku': (0.01, 1000.0), 'K': (1e-3, 1),
              'coop': (0.05, 20), 'Kdim': (0.5,2), 'omega': (0.05/20, 20/0.05)}
    return params

class DDDS(SG):
    def __init__(self, params={}):
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
            self.f3 = 0
            self.f12 = 0
            self.f23 = 0
            self.f13 = 0
            self.f123 = 0
            self.kb1 = 0
            self.ku1 = 0
            self.kb2 = 0
            self.ku2 = 0
            self.kb3 = 0
            self.ku3 = 0
            self.bcoop12 = 0
            self.ucoop12 = 0
            self.bcoop13 = 0
            self.ucoop13 = 0
            self.bcoop23 = 0
            self.ucoop23 = 0
            self.bcoop123 = 0
            self.ucoop123 = 0
            self.vol = 0
        else:
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
            self.f1 = params['f1'];
            self.f2 = params['f2'];
            self.f3 = params['f3'];
            self.f12 = params['f12'];
            self.f13 = params['f13'];
            self.f23 = params['f23'];
            self.f123 = params['f123'];
            self.kb1 = params['kb1'];
            self.kb2 = params['kb2'];
            self.kb3 = params['kb3'];
            if 'ku1' in params:
                self.ku1 = params['ku1'];
            elif 'Kd1' in params:
                self.ku1 = self.kb1 / params['Kd1'];
            if 'ku2' in params:
                self.ku2 = params['ku2'];
            elif 'Kd2' in params:
                self.ku2 = self.kb2 / params['Kd2'];
            if 'ku3' in params:
                self.ku3 = params['ku3'];
            elif 'Kd3' in params:
                self.ku3 = self.kb3 / params['Kd3'];

            self.bcoop12 = params['bcoop12'];
            self.bcoop13 = params['bcoop13'];
            self.bcoop23 = params['bcoop23'];
            self.bcoop123 = params['bcoop123'];
            if 'ucoop12' in params:
                self.ucoop12 = params['ucoop12'];
            elif 'omega12' in params:
                self.ucoop12 = params['omega12']/self.bcoop12;
            if 'ucoop23' in params:
                self.ucoop23 = params['ucoop23'];
            elif 'omega23' in params:
                self.ucoop23 = params['omega23']/self.bcoop23;
            if 'ucoop13' in params:
                self.ucoop13 = params['ucoop13'];
            elif 'omega13' in params:
                self.ucoop13 = params['omega13']/self.bcoop13;
            if 'ucoop123' in params:
                self.ucoop123 = params['ucoop123'];
            elif 'omega123' in params:
                self.ucoop123 = params['omega123']/self.bcoop123;
            self.vol = params['vol'];
        self.DNAstates = ['DNA0', 'DNA1', 'DNA2', 'DNA3', 'DNA12', 'DNA23', 'DNA13', 'DNA123'];
        self.allvars = np.array(['DNA0', 'DNA1', 'DNA2', 'DNA3', 'DNA12', 'DNA13', 'DNA23', 'mRNA', 'm', 'd']);
        self.nameparameters = np.array(['beta', 'gammam', 'gammamRNA', 'gammad', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3','f12', 'f13','f23','f123','kb1','kb2','kb3','ku1', 'ku2', 'ku3','bcoop12','bcoop13','bcoop23','bcoop123', 'ucoop12','ucoop13','ucoop23','ucoop123']);

        self.ss = {};  # steady state

    def outRange(self):
        physrange = physiologicalRange();
        line = "";
        par = [self.beta, self.gammam, self.gammamRNA, self.gammad, self.beta, self.alphaass, self.alphadiss, self.phi0,
               self.f1, self.f2, self.f3, self.f12, self.f13, self.f23, self.f123,
               self.kb1, self.ku1, self.kb2, self.ku2, self.kb3, self.ku3, self.bcoop12, self.ucoop12,
               self.bcoop13, self.ucoop13, self.bcoop23, self.ucoop23, self.bcoop123, self.ucoop123];
        name = ['beta', 'gammam', 'gammamRNA', 'gammad', 'beta', 'alphaass', 'alphadiss', 'phi0', 'f1', 'f2', 'f3', 'f12', 'f13', 'f23', 'f123',
                'kb1', 'ku1', 'kb2', 'ku2', 'kb3', 'ku3', 'bcoop12', 'ucoop12', 'bcoop13', 'ucoop13', 'bcoop23', 'ucoop23', 'bcoop123', 'ucoop123']
        parname = ['beta', 'gammam', 'gammamRNA', 'gammad', 'beta', 'alphaass', 'alphadiss', 'phi0', 'fa', 'fa', 'fa', 'fa', 'fa', 'fa', 'fr',
                   'kb', 'ku', 'kb', 'ku', 'kb', 'ku', 'bcoop', 'ucoop', 'bcoop', 'ucoop', 'bcoop', 'ucoop', 'bcoop', 'ucoop']
        for i in range(len(par)):
            if (par[i] < physrange[parname[i]][0]):
                line += "%s low; " % name[i];
            elif (par[i] == physrange[parname[i]][0]):
                line += "%s low bound; " % name[i];
            elif (par[i] > physrange[parname[i]][1]):
                line += "%s high; " % name[i];
            elif (par[i] == physrange[parname[i]][1]):
                line += "%s high bound; " % name[i];
        return line;

    def physicalizeSsLrpB(self):
        parname = ['beta', 'gammam', 'gammamRNA', 'gammad', 'beta', 'alphaass', 'alphadiss', 'phi0', 'f', 'f', 'f',
                   'f', 'f', 'f', 'f', 'kb', 'ku', 'kb', 'ku', 'kb', 'ku', 'bcoop', 'ucoop', 'bcoop', 'ucoop',
                   'bcoop', 'ucoop', 'bcoop', 'ucoop']
        physrange = physiologicalRange()
        ct = 1e-3 / 2.4

        Kd1 = 73.5 * ct
        Kd2 = 0.7 * ct
        Kd3 = 49.1 * ct
        omega12 = 4.1
        omega13 = 7.9
        omega23 = 2.1
        omega123 = 3.1
        if self.kb1 < physrange['kb'][0]:
            self.kb1 = physrange['kb'][0]
            self.ku1 = self.kb1/Kd1
        if self.kb1 > physrange['kb'][1]:
            self.kb1 = physrange['kb'][1]
            self.ku1 = self.kb1 / Kd1

        if self.ku1 < physrange['ku'][0]:
            self.ku1 = physrange['ku'][0]
            self.kb1 = self.ku1*Kd1
        if self.ku1 > physrange['ku'][1]:
            self.ku1 = physrange['ku'][1]
            self.kb1 = self.ku1 * Kd1

        if self.kb2 < physrange['kb'][0]:
            self.kb2 = physrange['kb'][0]
            self.ku2 = self.kb2 / Kd2
        if self.kb2 > physrange['kb'][1]:
            self.kb2 = physrange['kb'][1]
            self.ku2 = self.kb2 / Kd2

        if self.ku2 < physrange['ku'][0]:
            self.ku2 = physrange['ku'][0]
            self.kb2 = self.ku2 * Kd2
        if self.ku2 > physrange['ku'][1]:
            self.ku2 = physrange['ku'][1]
            self.kb2 = self.ku2 * Kd2

        if self.kb3 < physrange['kb'][0]:
            self.kb3 = physrange['kb'][0]
            self.ku3 = self.kb3 / Kd3
        if self.kb3 > physrange['kb'][1]:
            self.kb3 = physrange['kb'][1]
            self.ku3 = self.kb3 / Kd3

        if self.ku3 < physrange['ku'][0]:
            self.ku3 = physrange['ku'][0]
            self.kb3 = self.ku3 * Kd3
        if self.ku3 > physrange['ku'][1]:
            self.ku3 = physrange['ku'][1]
            self.kb3 = self.ku3 * Kd3

        if self.bcoop12 < physrange['coop'][0]:
            self.bcoop12 = physrange['coop'][0]
            self.ucoop12 = self.bcoop12 / omega12
        if self.bcoop12 > physrange['coop'][1]:
            self.bcoop12 = physrange['coop'][1]
            self.ucoop12 = self.bcoop12 / omega12

        if self.ucoop12 < physrange['coop'][0]:
            self.ucoop12 = physrange['coop'][0]
            self.bcoop12 = self.ucoop12 * omega12
        if self.ucoop12 > physrange['coop'][1]:
            self.ucoop12 = physrange['coop'][1]
            self.bcoop12 = self.ucoop12 * omega12

        if self.bcoop13 < physrange['coop'][0]:
            self.bcoop13 = physrange['coop'][0]
            self.ucoop13 = self.bcoop13 / omega13
        if self.bcoop13 > physrange['coop'][1]:
            self.bcoop13 = physrange['coop'][1]
            self.ucoop13 = self.bcoop13 / omega13

        if self.ucoop13 < physrange['coop'][0]:
            self.ucoop13 = physrange['coop'][0]
            self.bcoop13 = self.ucoop13 * omega13
        if self.ucoop13 > physrange['coop'][1]:
            self.ucoop13 = physrange['coop'][1]
            self.bcoop13 = self.ucoop13 * omega13

        if self.bcoop23 < physrange['coop'][0]:
            self.bcoop23 = physrange['coop'][0]
            self.ucoop23 = self.bcoop23 / omega23
        if self.bcoop23 > physrange['coop'][1]:
            self.bcoop23 = physrange['coop'][1]
            self.ucoop23 = self.bcoop23 / omega23

        if self.ucoop23 < physrange['coop'][0]:
            self.ucoop23 = physrange['coop'][0]
            self.bcoop23 = self.ucoop23 * omega23
        if self.ucoop23 > physrange['coop'][1]:
            self.ucoop23 = physrange['coop'][1]
            self.bcoop23 = self.ucoop23 * omega23

        if self.bcoop123 < physrange['coop'][0]:
            self.bcoop123 = physrange['coop'][0]
            self.ucoop123 = self.bcoop123 / omega123
        if self.bcoop123 > physrange['coop'][1]:
            self.bcoop123 = physrange['coop'][1]
            self.ucoop123 = self.bcoop123 / omega123

        if self.ucoop123 < physrange['coop'][0]:
            self.ucoop123 = physrange['coop'][0]
            self.bcoop123 = self.ucoop123 * omega123
        if self.ucoop123 > physrange['coop'][1]:
            self.ucoop123 = physrange['coop'][1]
            self.bcoop123 = self.ucoop123 * omega123

    # return string of paramters
    def parameterline(self):
        string = "%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n" % (
            self.beta, self.gammam, self.gammamRNA, self.gammad, self.alphaass, self.alphadiss, self.phi0,
            self.f1, self.f2, self.f3, self.f12, self.f13, self.f23, self.f123,
            self.kb1, self.kb2, self.kb3, self.ku1, self.ku2, self.ku3, self.bcoop12, self.bcoop13,
            self.bcoop23, self.bcoop123, self.ucoop12, self.ucoop13, self.ucoop23, self.ucoop123, self.taum,
            self.taumRNA);
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
        self.f1 = params['f1'];
        self.f2 = params['f2'];
        self.f3 = params['f3'];
        self.f12 = params['f12'];
        self.f13 = params['f13'];
        self.f23 = params['f23'];
        self.f123 = params['f123'];
        self.kb1 = params['kb1'];
        self.kb2 = params['kb2'];
        self.kb3 = params['kb3'];
        if 'ku1' in params:
            self.ku1 = params['ku1'];
            self.Kd1 = self.kb1/self.ku1
        elif 'Kd1' in params:
            self.Kd1 = params['Kd1'];
            self.ku1 = self.kb1/self.Kd1;
        if 'ku2' in params:
            self.ku2 = params['ku2'];
            self.Kd2 = self.kb2/self.ku2
        elif 'Kd2' in params:
            self.Kd2 = params['Kd2'];
            self.ku2 = self.kb2 / self.Kd2;
        if 'ku3' in params:
            self.ku3 = params['ku3'];
            self.Kd3 = self.kb3/self.ku3
        elif 'Kd3' in params:
            self.Kd3 = params['Kd3'];
            self.ku3 = self.kb3 / self.Kd3;

        self.bcoop12 = params['bcoop12'];
        self.bcoop13 = params['bcoop13'];
        self.bcoop23 = params['bcoop23'];
        self.bcoop123 = params['bcoop123'];
        self.ucoop12 = params['ucoop12'];
        self.ucoop13 = params['ucoop13'];
        self.ucoop23 = params['ucoop23'];
        self.ucoop123 = params['ucoop123'];
        self.vol = params['vol'];
        self.ss.clear();

    # set parameters to mutant of parent
    def smallmutant(self, parent, constraints, withdelay=False, SsLrpB=True):
        if (withdelay):
            idx = int(rd.uniform(0, 30));
        else:
            idx = int(rd.uniform(0, 28));
        mutation = np.ones(30);
        mutation[idx] = 10 ** (rd.uniform(np.log10(0.5), np.log10(2)));
        if(len(constraints)==0):
            print('no constraints')
        else:
            self.beta = max(min(parent.beta * mutation[0], constraints['beta'][1]), constraints['beta'][0]);
            self.gammam = max(min(parent.gammam * mutation[1], constraints['gammam'][1]), constraints['gammam'][0]);
            self.gammamRNA = max(min(parent.gammamRNA * mutation[2], constraints['gammamRNA'][1]), constraints['gammamRNA'][0]);
            self.gammad = max(min(parent.gammad * mutation[3], constraints['gammad'][1]), constraints['gammad'][0]);
            self.alphaass = max(min(parent.alphaass * mutation[4], constraints['alphaass'][1]), constraints['alphaass'][0]);
            Kdim = max(min(self.alphadiss/self.alphaass * mutation[5], constraints['Kdim'][1]), constraints['Kdim'][0]);
            self.alphadiss = self.alphaass*Kdim;
            self.phi0 = max(min(parent.phi0 * mutation[6], constraints['phi0'][1]), constraints['phi0'][0]);
            self.f1 = max(min(parent.f1 * mutation[7], constraints['f'][1]), constraints['f'][0]); #todo fa
            self.f2 = max(min(parent.f2 * mutation[8], constraints['f'][1]), constraints['f'][0]);#todo fa
            self.f3 = max(min(parent.f3 * mutation[9], constraints['f'][1]), constraints['f'][0]);#todo fa
            self.f12 = max(min(parent.f12 * mutation[10], constraints['f'][1]), constraints['f'][0]);#todo fa
            self.f13 = max(min(parent.f13 * mutation[11], constraints['f'][1]), constraints['f'][0]);#todo fa
            self.f23 = max(min(parent.f23 * mutation[12], constraints['f'][1]), constraints['f'][0]);#todo fa
            #self.f123 = max(min(parent.f123 * mutation[13], 0.5*max(self.f1, self.f2, self.f3, self.f12, self.f23, self.f23)), constraints['f'][0]);#todo fr
            self.f123 = max(min(parent.f123 * mutation[13], constraints['f'][1]), constraints['f'][0]);#todo fr
            self.kb1 = max(min(parent.kb1 * mutation[14], constraints['kb'][1]), constraints['kb'][0]);
            self.ku1 = max(min(parent.ku1 * mutation[15], constraints['ku'][1]), constraints['ku'][0]);
            self.kb2 = max(min(parent.kb2 * mutation[16], constraints['kb'][1]), constraints['kb'][0]);
            self.ku2 = max(min(parent.ku2 * mutation[17], constraints['ku'][1]), constraints['ku'][0]);
            self.kb3 = max(min(parent.kb3 * mutation[18], constraints['kb'][1]), constraints['kb'][0]);
            self.ku3 = max(min(parent.ku3 * mutation[19], constraints['ku'][1]), constraints['ku'][0]);
            self.bcoop12 = max(min(parent.bcoop12 * mutation[20], constraints['bcoop'][1]), constraints['bcoop'][0]);
            self.ucoop12 = max(min(parent.ucoop12 * mutation[21], constraints['ucoop'][1]), constraints['ucoop'][0]);
            self.bcoop23 = max(min(parent.bcoop23 * mutation[22], constraints['bcoop'][1]), constraints['bcoop'][0]);
            self.ucoop23 = max(min(parent.ucoop23 * mutation[23], constraints['ucoop'][1]), constraints['ucoop'][0]);
            self.bcoop13 = max(min(parent.bcoop13 * mutation[24], constraints['bcoop'][1]), constraints['bcoop'][0]);
            self.ucoop13 = max(min(parent.ucoop123 * mutation[25], constraints['ucoop'][1]), constraints['ucoop'][0]);
            self.bcoop123 = max(min(parent.bcoop123 * mutation[26], constraints['bcoop'][1]), constraints['bcoop'][0]);
            self.ucoop123 = max(min(parent.ucoop123 * mutation[27], constraints['ucoop'][1]), constraints['ucoop'][0]);
            self.taum = min(5, parent.taum * mutation[28]);
            self.taumRNA = min(5, parent.taumRNA * mutation[29]);
            self.vol = parent.vol;
        if(SsLrpB):
            self.ku1 = self.kb1 / 1.764
            self.ku2 = self.kb2 / 0.0168
            self.ku3 = self.kb3 / 1.1784
            self.ucoop12 = self.bcoop12 / 4.1
            self.ucoop23 = self.bcoop23 / 7.9
            self.ucoop13 = self.bcoop13 / 2.1
            self.ucoop123 = self.bcoop123 / 3.1
        self.ss.clear();
        return

    def largemutant(self, parent, constraints, withdelay=False, SsLrpB=True):
        idx = int(rd.uniform(0, 28))
        self.copy(parent)
        if (idx == 0):
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
            self.f1 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo f
            #self.f123 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(0.5*max(self.f1, self.f2, self.f3, self.f12, self.f13, self.f23)))) #todo f
        elif (idx == 8):
            self.f2 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo f
            #self.f123 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(0.5*max(self.f1, self.f2, self.f3, self.f12, self.f13, self.f23)))) #todo f
        elif (idx == 9):
            self.f3 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo f
            #self.f123 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(0.5*max(self.f1, self.f2, self.f3, self.f12, self.f13, self.f23)))) #todo f
        elif (idx == 10):
            self.f12 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo f
            #self.f123 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(0.5*max(self.f1, self.f2, self.f3, self.f12, self.f13, self.f23)))) #todo f
        elif (idx == 11):
            self.f13 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo f
            #self.f123 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(0.5*max(self.f1, self.f2, self.f3, self.f12, self.f13, self.f23)))) #todo f
        elif (idx == 12):
            self.f23 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo f
            #self.f123 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(0.5*max(self.f1, self.f2, self.f3, self.f12, self.f13, self.f23)))) #todo f
        elif (idx == 13):
            self.f123 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(constraints['f'][1]))) #todo f
            #self.f123 = np.exp(rd.uniform(np.log(constraints['f'][0]), np.log(0.5*max(self.f1, self.f2, self.f3, self.f12, self.f13, self.f23)))) #todo f
        elif (idx == 14):
            self.kb1 = np.exp(rd.uniform(np.log(constraints['kb'][0]), np.log(constraints['kb'][1])))
        elif (idx == 15):
            self.kb2 = np.exp(rd.uniform(np.log(constraints['kb'][0]), np.log(constraints['kb'][1])))
        elif (idx == 16):
            self.kb3 = np.exp(rd.uniform(np.log(constraints['kb'][0]), np.log(constraints['kb'][1])))
        elif (idx == 17):
            self.ku1 = np.exp(rd.uniform(np.log(constraints['ku'][0]), np.log(constraints['ku'][1])))
        elif (idx == 18):
            self.ku2 = np.exp(rd.uniform(np.log(constraints['ku'][0]), np.log(constraints['ku'][1])))
        elif (idx == 19):
            self.ku3 = np.exp(rd.uniform(np.log(constraints['ku'][0]), np.log(constraints['ku'][1])))
        elif (idx == 20):
            self.bcoop12 = np.exp(rd.uniform(np.log(constraints['bcoop'][0]), np.log(constraints['bcoop'][1])))
        elif (idx == 21):
            self.bcoop13 = np.exp(rd.uniform(np.log(constraints['bcoop'][0]), np.log(constraints['bcoop'][1])))
        elif (idx == 22):
            self.bcoop23 = np.exp(rd.uniform(np.log(constraints['bcoop'][0]), np.log(constraints['bcoop'][1])))
        elif (idx == 23):
            self.bcoop123 = np.exp(rd.uniform(np.log(constraints['bcoop'][0]), np.log(constraints['bcoop'][1])))
        elif (idx == 24):
            self.ucoop12 = np.exp(rd.uniform(np.log(constraints['ucoop'][0]), np.log(constraints['ucoop'][1])))
        elif (idx == 25):
            self.ucoop13 = np.exp(rd.uniform(np.log(constraints['ucoop'][0]), np.log(constraints['ucoop'][1])))
        elif (idx == 26):
            self.ucoop23 = np.exp(rd.uniform(np.log(constraints['ucoop'][0]), np.log(constraints['ucoop'][1])))
        elif (idx == 27):
            self.ucoop123 = np.exp(rd.uniform(np.log(constraints['ucoop'][0]), np.log(constraints['ucoop'][1])))
        self.ss.clear();

    def copy(self, parent):
        self.beta = parent.beta
        self.gammam = parent.gammam
        self.gammamRNA = parent.gammamRNA
        self.gammad = parent.gammad
        self.alphaass = parent.alphaass
        self.alphadiss = parent.alphadiss
        self.DNAtot = parent.DNAtot;
        self.phi0 = parent.phi0
        self.f1 = parent.f1
        self.f2 = parent.f2
        self.f3 = parent.f3
        self.f12 = parent.f12
        self.f13 = parent.f13
        self.f23 = parent.f23
        self.f123 = parent.f123
        self.kb1 = parent.kb1
        self.kb2 = parent.kb2
        self.kb3 = parent.kb3
        self.ku1 = parent.ku1
        self.ku2 = parent.ku2
        self.ku3 = parent.ku3
        self.bcoop12 = parent.bcoop12
        self.bcoop13 = parent.bcoop13
        self.bcoop23 = parent.bcoop23
        self.bcoop123 = parent.bcoop123
        self.ucoop12 = parent.ucoop12
        self.ucoop13 = parent.ucoop13
        self.ucoop23 = parent.ucoop23
        self.ucoop123 = parent.ucoop123
        self.taum = parent.taum
        self.taumRNA = parent.taumRNA
        self.vol = parent.vol;

    def setsteadystate(self):
        fdenom = lambda d: (self.kb1 * self.kb2 * self.kb3 / (self.ku1 * self.ku2 * self.ku3) * self.bcoop12 * self.bcoop23 * self.bcoop13 * self.bcoop123 / (self.ucoop12 * self.ucoop23 * self.ucoop13 * self.ucoop123) * d ** 3 + (self.kb1 * self.kb2 / (self.ku1 * self.ku2) * self.bcoop12 / self.ucoop12 + self.kb1 * self.kb3 / (self.ku1 * self.ku3) * self.bcoop13 / self.ucoop13 + self.kb2 * self.kb3 / (self.ku2 * self.ku3) * self.bcoop23 / self.ucoop23) * d ** 2 + (self.kb1 / self.ku1 + self.kb2 / self.ku2 + self.kb3 / self.ku3) * d + 1)
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

        func = lambda x: (self.phi0 * ((1 - self.f123)*fDNA0(fd(x)) + (self.f1 - self.f123) * fDNA1(fd(x))
                                       + (self.f2 - self.f123)*fDNA2(fd(x)) + (self.f3 - self.f123)*fDNA3(fd(x))
                        + (self.f12 - self.f123)*fDNA12(fd(x)) + (self.f13 - self.f123)*fDNA13(fd(x))
                        + (self.f23 - self.f123)*fDNA23(fd(x)) + self.f123 * self.DNAtot)
                          - self.gammamRNA / self.beta * (
                              self.gammam * x) - self.gammamRNA / self.beta * 2 * self.alphaass * x ** 2 / (
                              self.alphadiss / self.gammad + 1))

        if (func(1e25) > 0):
            m = np.array([np.inf])
        elif (func(1e-25) < 0):
            m = np.array([0])
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

        if (False):
            x = np.logspace(-19, 25, 500)
            plt.figure()
            plt.plot(x, func(x))
            plt.xscale('log')
            plt.ylim([-0.1, 0.1])
            plt.plot(m, 0, '*')
            plt.show()

        d = fd(m);
        DNA0 = fDNA0(d);
        DNA1 = fDNA1(d);
        DNA2 = fDNA2(d);
        DNA3 = fDNA3(d);
        DNA12 = fDNA12(d);
        DNA13 = fDNA13(d);
        DNA23 = fDNA23(d);
        DNA123 = self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23;
        mRNA = self.phi0 * ((1 - self.f123) * DNA0 + (self.f1 - self.f123) * DNA1
                            + (self.f2 - self.f123) * DNA2 + (self.f3 - self.f123) * DNA3
                            + (self.f12 - self.f123) * DNA12 + (self.f13 - self.f123) * DNA13
                            + (self.f23 - self.f123) * DNA23 + self.f123 * self.DNAtot) / self.gammamRNA

        self.ss = {'DNA0': DNA0, 'DNA1': DNA1, 'DNA2': DNA2, 'DNA3': DNA3, 'DNA12': DNA12, 'DNA13': DNA13,
                   'DNA23': DNA23, 'DNA123': DNA123, 'mRNA': mRNA, 'm': m, 'd': d}

    def eqns(self, var, t, ts, qs):
        DNA0, DNA1, DNA2, DNA3, DNA12, DNA13, DNA23, mRNA, m, d = var;
        # ts, qs = args;

        dDNA0 = (self.ku1 * DNA1 + self.ku2 * DNA2 + self.ku3 * DNA3) - (
                                                                               self.kb1 + self.kb2 + self.kb3) * d * DNA0
        dDNA1 = self.kb1 * d * DNA0 + self.ucoop12 * self.ku2 * DNA12 + self.ucoop13 * self.ku3 * DNA13 - (
                                                                                                                 self.ku1 + (
                                                                                                                     self.bcoop12 * self.kb2 + self.bcoop13 * self.kb3) * d) * DNA1
        dDNA2 = self.kb2 * d * DNA0 + self.ucoop12 * self.ku1 * DNA12 + self.ucoop23 * self.ku3 * DNA23 - (
                                                                                                                 self.ku2 + (
                                                                                                                     self.bcoop12 * self.kb1 + self.bcoop23 * self.kb3) * d) * DNA2
        dDNA3 = self.kb3 * d * DNA0 + self.ucoop13 * self.ku1 * DNA13 + self.ucoop23 * self.ku2 * DNA23 - (
                                                                                                                 self.ku3 + (
                                                                                                                     self.bcoop13 * self.kb1 + self.bcoop23 * self.kb2) * d) * DNA3
        dDNA12 = self.bcoop12 * d * (
            self.kb1 * DNA2 + self.kb2 * DNA1) + self.ucoop123 * self.ucoop13 * self.ucoop23 * self.ku3 * (
            self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (
                                                                                   self.bcoop123 * self.bcoop13 * self.bcoop23 * self.kb3 * d + self.ucoop12 * (
                                                                                       self.ku1 + self.ku2)) * DNA12
        dDNA13 = self.bcoop13 * d * (
            self.kb1 * DNA3 + self.kb3 * DNA1) + self.ucoop123 * self.ucoop12 * self.ucoop23 * self.ku2 * (
            self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (
                                                                                   self.bcoop123 * self.bcoop12 * self.bcoop23 * self.kb2 * d + self.ucoop13 * (
                                                                                       self.ku1 + self.ku3)) * DNA13
        dDNA23 = self.bcoop23 * d * (
            self.kb3 * DNA2 + self.kb2 * DNA3) + self.ucoop123 * self.ucoop12 * self.ucoop13 * self.ku1 * (
            self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13) - (
                                                                                   self.bcoop123 * self.bcoop12 * self.bcoop13 * self.kb1 * d + self.ucoop23 * (
                                                                                       self.ku2 + self.ku3)) * DNA23
        dmRNA = self.phi0 * (DNA0 + self.f1 * DNA1 + self.f2 * DNA2 + self.f3 * DNA3
                             + self.f12 * DNA12 + self.f23 * DNA23 + self.f13 *DNA13 + self.f123 * (
            self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23)) - self.gammamRNA * mRNA
        dm = self.beta * mRNA - 2 * self.alphaass * pow(m, 2) + 2 * self.alphadiss * d - self.gammam * m
        dd = (-d * (
        (self.kb1 + self.kb2 + self.kb3) * DNA0 + (self.kb2 * self.bcoop12 + self.kb3 * self.bcoop13) * DNA1 + (
            self.kb1 * self.bcoop12 + self.kb3 * self.bcoop23) * DNA2 + (
            self.kb1 * self.bcoop13 + self.kb2 * self.bcoop23) * DNA3 + self.kb3 * self.bcoop123 * self.bcoop13 * self.bcoop23 * DNA12 + self.kb2 * self.bcoop123 * self.bcoop12 * self.bcoop23 * DNA13 + self.kb1 * self.bcoop123 * self.bcoop12 * self.bcoop13 * DNA23) + (
                  self.ku1 * DNA1 + self.ku2 * DNA2 + self.ku3 * DNA3 + (
                      self.ku1 + self.ku2) * self.ucoop12 * DNA12 + (self.ku1 + self.ku3) * self.ucoop13 * DNA13 + (
                      self.ku2 + self.ku3) * self.ucoop23 * DNA23 + self.ucoop123 * (
                      self.ku1 * self.ucoop12 * self.ucoop13 + self.ku2 * self.ucoop12 * self.ucoop23 + self.ku3 * self.ucoop13 * self.ucoop23) * (
                      self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA13 - DNA23)) + self.alphaass * pow(m,
                                                                                                              2) - self.alphadiss * d - self.gammad * d)
        return [dDNA0, dDNA1, dDNA2, dDNA3, dDNA12, dDNA13, dDNA23, dmRNA, dm, dd];

    # delay differential equations in dictionary format for pydelay
    def eqnsdelay(self):
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

    def histfunc(self, qs):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        hf = {'m': lambda t: self.ss['m'] / qs * 1.2, 'd': lambda t: self.ss['d'] / qs * 1.2,
              'mRNA': lambda t: self.ss['mRNA'] / qs * 1.2, 'DNA0': lambda t: self.ss['DNA0'] / qs,
              'DNA1': lambda t: self.ss['DNA1'] / qs, 'DNA2': lambda t: self.ss['DNA2'] / qs,
              'DNA3': lambda t: self.ss['DNA3'] / qs, 'DNA12': lambda t: self.ss['DNA12'] / qs,
              'DNA13': lambda t: self.ss['DNA13'] / qs, 'DNA23': lambda t: self.ss['DNA23'] / qs, }
        return hf

    def quasisteadystate(self, x, monomer=True):
        if (monomer):
            m = x;
            d = self.alphaass * m ** 2 / (self.alphadiss + self.gammad)
        else:
            d = x;
            m = np.sqrt((self.alphadiss + self.gammad) * d / self.alphaass);
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
        ss = {'DNA0': DNA0, 'DNA1': DNA1, 'DNA2': DNA2, 'DNA3': DNA3, 'DNA12': DNA12, 'DNA13': DNA13, 'DNA23': DNA23,
              'DNA123': DNA123, 'd': d, 'm': m}
        return ss

    def transcriptionrate(self, x, monomer=True):
        ssDNA = self.quasisteadystate(x, monomer)

        return self.phi0 * (ssDNA['DNA0'] + self.f1 * ssDNA['DNA1'] + self.f2 * ssDNA['DNA2']
                            + self.f3 * ssDNA['DNA3'] + self.f12  * ssDNA['DNA12']
                            + self.f23 * ssDNA['DNA23'] + self.f13  * ssDNA['DNA13'] + self.f123 * (ssDNA['DNA123']))

    def translationrate(self, x, monomer=True):
        return self.beta*self.transcriptionrate(x, monomer)/self.gammamRNA

    def nullClineM(self, x, monomer=True):  # TODO: IF NOT MONOMER RETURN DIMER STATE???
        qss = self.quasisteadystate(x, monomer);
        # values of mRNA for which dm/dt are 0 :
        m0 = (2 * self.alphaass * pow(qss['m'], 2) - 2 * self.alphadiss * qss['d'] + self.gammam * qss['m']) / self.beta
        return m0;

    def nullClineMRNA(self, x, monomer=True):
        qss = self.quasisteadystate(x, monomer);
        # values of mRNA for which dmRNA/dt are 0 :
        mRNA0 = self.phi0 * ((1 - self.f123) * qss['DNA0'] + (self.f1 - self.f123) * qss['DNA1']
                             + (self.f2 - self.f123) * qss['DNA2'] + (self.f3 - self.f123) * qss['DNA3']
                             + (self.f12 - self.f123) * qss['DNA12'] + (self.f13 - self.f123) * qss['DNA13']
                             + (self.f23 - self.f123) * qss['DNA23'] + self.f123 * self.DNAtot) / self.gammamRNA
        return mRNA0;

    def paramsnd(self, ts, qs):
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

    def timelapse(self, t_f, dt, t_i=0, hf={}):
        if (not len(self.ss) > 0):
            self.setsteadystate();
        ts = 1.0  # 1.0/self.gammam;
        qs = 1.0  # self.ss['m']
        if (self.taum == 0 and self.taumRNA == 0):
            if (len(hf) > 0):
                y0 = [hf['DNA0'](0), hf['DNA1'](0), hf['DNA2'](0), hf['DNA3'](0), hf['DNA12'](0), hf['DNA13'](0),
                      hf['DNA23'](0), hf['mRNA'](0), hf['m'](0), hf['d'](0)]
            else:
                #y0 = [min(self.ss['DNA0']), min(self.ss['DNA1']), min(self.ss['DNA2']), min(self.ss['DNA3']), min(self.ss['DNA12']),
                #      min(self.ss['DNA13']), min(self.ss['DNA23']), min(self.ss['mRNA']) * 1.2, min(self.ss['m']) * 1.2, min(self.ss['d']) * 1.2]
                y0 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            t = np.arange(t_i, t_f, dt)
            params = {'qs': qs, 'ts': ts}
            ode = ode23(self, odeint(self.eqns, y0, t, args=(ts, qs)), t, []);
            return ode
        else:
            dde = dde23(eqns=self.eqnsdelay(), params=self.paramsnd(ts, qs))
            dde.set_sim_params(tfinal=t_f, dtmax=dt, AbsTol=10 ** -6, RelTol=10 ** -3)
            if (len(hf) > 0):
                dde.hist_from_funcs(hf, 50)
            else:
                dde.hist_from_funcs(self.histfunc(qs), 50)
            dde.run()
            return dde

    def Jacobian(self, var):
        m, d, DNA0, DNA1, DNA2, DNA3, DNA12, DNA13, DNA23 = var;
        c12DNA123 = self.ucoop123 * self.ucoop13 * self.ucoop23 * self.ku3
        c13DNA123 = self.ucoop123 * self.ucoop12 * self.ucoop23 * self.ku2
        c23DNA123 = self.ucoop123 * self.ucoop12 * self.ucoop13 * self.ku1
        cmDNA123 = self.ucoop123 * (
            self.ku1 * self.ucoop12 * self.ucoop13 + self.ku2 * self.ucoop12 * self.ucoop23 + self.ku3 * self.ucoop13 * self.ucoop23)

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

    def reactionmatrix(self):
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
        DNA0, DNA1, DNA2, DNA3, DNA12, DNA23, DNA13, mRNA, m, d = var;
        DNA123 = self.DNAtot - DNA0 - DNA1 - DNA2 - DNA3 - DNA12 - DNA23 - DNA13;
        prop = np.array(
                  [self.kb1*DNA0*d,
                    self.ku1*DNA1,
                    self.kb2*DNA0*d,
                    self.ku2*DNA2,
                    self.kb3 * DNA0 * d,
                    self.ku3 * DNA3,
                    self.kb2*self.bcoop12*DNA1*d,
                    self.ku2*self.ucoop12*DNA12,
                    self.kb3 * self.bcoop13 * DNA1 * d,
                    self.ku3 * self.ucoop13 * DNA13,
                    self.kb1 * self.bcoop12 * DNA2 * d,
                    self.ku1 * self.ucoop12 * DNA12,
                    self.kb3 * self.bcoop23 * DNA2 * d,
                    self.ku3 * self.ucoop23 * DNA23,
                    self.kb1 * self.bcoop13 * DNA3 * d,
                    self.ku1 * self.ucoop13 * DNA13,
                    self.kb2 * self.bcoop23 * DNA3 * d,
                    self.ku2 * self.ucoop23 * DNA23,
                    self.kb3 * self.bcoop13 * self.bcoop23 * self.bcoop123 * DNA12*d,
                    self.ku3 * self.ucoop13 * self.ucoop23 * self.ucoop123 * DNA123,
                    self.kb1 * self.bcoop13 * self.bcoop12 * self.bcoop123 * DNA23 * d,
                    self.ku1 * self.ucoop13 * self.ucoop12 * self.ucoop123 * DNA123,
                    self.kb2 * self.bcoop12 * self.bcoop23 * self.bcoop123 * DNA13 * d,
                    self.ku2 * self.ucoop12 * self.ucoop23 * self.ucoop123 * DNA123,
                    self.phi0*(DNA0 + self.f1*DNA1 + self.f2*DNA2 + self.f3*DNA3 + self.f12*DNA12 + self.f13*DNA13 + self.f23*DNA23 + self.f123*DNA123),
                    self.beta*mRNA,
                    self.alphaass*m*(m-1),
                    self.alphadiss*d,
                    self.gammamRNA*mRNA,
                    self.gammam*m,
                    self.gammad*d])
        if(np.any(prop<0)):
            print(prop)
        return prop

    #TODO return right value
    def eigenvalues(self, ax=0, color='k'):
        if (not len(self.ss) > 0):
            self.steadystate();
        if (max(self.ss['m']) > 1e30):
            return [np.nan]

        for i in range(len(self.ss['m'])):
            var = [self.ss['m'][i], self.ss['d'][i], self.ss['DNA0'][i], self.ss['DNA1'][i], self.ss['DNA2'][i], self.ss['DNA3'][i],
                   self.ss['DNA12'][i], self.ss['DNA13'][i], self.ss['DNA23'][i]];
            J = self.Jacobian(var);
            eigvals = np.linalg.eigvals(J)
            if (ax != 0):
                for eigval in eigvals:
                    ax.plot(eigval.real, eigval.imag, c=color, marker='o')
                ax.axvline(0, color='k', linestyle=':')
                ax.set_xlabel("Real")
                ax.set_ylabel("Imaginary")
                ax.set_yticks([0])
        return eigvals