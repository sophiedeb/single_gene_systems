#!/usr/bin/env python

"""find_hull_3DS.py: Function to find hull of physiological range for 3DS."""

__author__ = "Lana Descheemaeker"
__email__ = "Lana.Descheemaeker@vub.be"

import multiprocessing
from functools import partial
import time
import numpy as np

num_cores = 4

def min_max_d(e=0, f=0):
    ''' Find the possible minimum and maximum of parameter D given E and F'''

    file = 'bistability/3DS/hull/%.3E.txt' % f

    # Set random seed
    s = (int(time.time())+int(1e4*f)+int(1e7*e))%157862
    np.random.seed(s)

    gooditerations = 0
    tries = 0

    # Initialization of the parameters
    minD = 1e23
    maxD = 0

    minKd1 = 0
    minKd2 = 0
    minKd3 = 0
    minomega12 = 0
    minomega13 = 0
    minomega23 = 0

    maxKd1 = 0
    maxKd2 = 0
    maxKd3 = 0
    maxomega12 = 0
    maxomega13 = 0
    maxomega23 = 0

    # Draw random numbers for the variables that dictate E and F, calculate D and store if the value is higher than
    # the maximum or lower than the minimum
    while (gooditerations < 1e3 and tries < 3e3):
        # Draw random Kds
        Kd1 = 10 ** np.random.uniform(np.log10(f / 3), min(4, np.log10(f-2e-4)))
        Kd2 = 10 ** np.random.uniform(-4, min(4, np.log10(f - Kd1-1e-4)))
        Kd3 = f - Kd1 - Kd2

        # Check if Kds are in physiological range and if good values to reconstruct E, else calculate new Kds
        if (1e4 > Kd3 > 1e-4 and 0.05/20*(Kd1*Kd2 + Kd2*Kd3 + Kd1*Kd3) < e < 20/0.05*(Kd1*Kd2 + Kd2*Kd3 + Kd1*Kd3)):
            # Order Kd values from high to low (Kd1 > Kd2 > Kd3)
            if (Kd1 < Kd2):
                temp = Kd1
                Kd1 = Kd2
                Kd2 = temp
            if (Kd2 < Kd3):
                temp = Kd2
                Kd2 = Kd3
                Kd3 = temp
            if (Kd1 < Kd2):
                temp = Kd1
                Kd1 = Kd2
                Kd2 = temp

            # Draw random omega12 & omega13
            omega12 = 10 ** np.random.uniform(np.log10(max(0.05/20, (e - 20/0.05*(Kd1*Kd3 + Kd2*Kd3))/(Kd1*Kd2))),
                                              np.log10(min(20/0.05, (e - 0.05/20*(Kd1*Kd3 + Kd2*Kd3))/(Kd1*Kd2))))
            omega13 = 10 ** np.random.uniform(np.log10(max(0.05/20, (e - (omega12*Kd1*Kd2 + 20/0.05*Kd2*Kd3))/(Kd1*Kd3))),
                                              np.log10(min(20/0.05, (e - (omega12*Kd1*Kd2 + 0.05/20*Kd2*Kd3))/(Kd1*Kd3))))
            omega23 = (e - omega12*Kd1*Kd2 - omega13*Kd1*Kd3)/(Kd2*Kd3)

            # Check if omega23 is in the physiological range otherwise start over by drawing new Kds
            if (0.05/20 <= omega23 <= 20/0.05):
                # Calculate D with the random parameters
                D = omega12*omega23*omega13*Kd1*Kd2*Kd3
                gooditerations += 1

                # If new D smaller than minimum or higher than maximum, store the value
                if D < minD:
                    minD = D
                    minKd1 = Kd1
                    minKd2 = Kd2
                    minKd3 = Kd3
                    minomega12 = omega12
                    minomega13 = omega13
                    minomega23 = omega23
                if D > maxD:
                    maxD = D
                    maxKd1 = Kd1
                    maxKd2 = Kd2
                    maxKd3 = Kd3
                    maxomega12 = omega12
                    maxomega13 = omega13
                    maxomega23 = omega23
            else:
                tries += 1
        else:
            tries += 1

    # Save the extrema
    data = (e, minD, minKd1, minKd2, minKd3, minomega12, minomega23, minomega13,
                         maxD, maxKd1, maxKd2, maxKd3, maxomega12, maxomega23, maxomega13)
    with open(file, 'a') as f_:
        f_.write('%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n'%data)

def main():
    p = multiprocessing.Pool(num_cores)

    # for every value of E and F find the minimum and maximum possible value of D
    for f in np.logspace(np.log10(3e-4), np.log10(3e4),NN):
        if(f-2e-4 <= 1e4):
            minx = 1e-8 + 2 * 1e-4 * (f - 2e-4)
        elif(f-1e-4-1e4 <= 1e4):
            minx = (1e-4*(f-1e4-1e-4) + 1e4*(f-1e4-1e-4) + 1)
        else:
            minx = ((f - 2e4)**2 + 1e4 * (f - 2e4) + 1e8)
        e = np.logspace(np.log10(0.05 / 20) + np.log10(minx),
                np.log10(f ** 2 / 3) + np.log10(20 / 0.05), NN)
        min_max_d_ = partial(min_max_d, f=f)
        p.map(min_max_d_, e)

if __name__ == "__main__":
    main()

