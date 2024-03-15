#!/usr/bin/env python3

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc

# TeX labels
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def pois(la, x):
    return np.exp(-la) * (la ** x) / math.gamma(x + 1)

numdeaths = np.asarray([1.0 * i for i in range(0, 5)])
numcases = np.asarray([109, 65, 22, 3, 1])
lam = 0.610
poispred = np.asarray([np.sum(numcases) * pois(lam, x) for x in numdeaths])

avgd = np.sum(numcases * numdeaths) / np.sum(numcases)

print("Num: ", numdeaths)
print("Data: ", numcases)
print("Poisson: ", poispred)
print("Average: ", avgd)

plt.bar(numdeaths, numcases, width = 0.5, color='b')
plt.bar(numdeaths + 0.5, poispred, width = 0.5, color='r')

plt.xlabel(r"deaths")
plt.ylabel(r"cases")

plt.savefig("ex8.1.eps", dpi=800)


#calculate chi2

def chi2(x, y):
    z = (x-y)**2 / y
    return np.sum(z)

print("Chi2 = ", chi2(numcases, poispred))
print("DoF = ", len(numcases) - 1)
