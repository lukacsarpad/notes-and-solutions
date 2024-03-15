#!/usr/bin/env python3

import numpy as np
import math
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib import rc

# TeX labels
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def pois(la, x):
    return np.exp(-la) * (la ** x) / math.gamma(x + 1)

numevs = np.asarray([1.0 * i for i in range(0, 10)])
numcases = np.asarray([1042, 860, 307, 78, 15, 3, 0, 0, 0, 1])

# without outliers
numcwooutl = np.copy(numcases)
numcwooutl[9] = 0

avg = np.sum(numevs * numcases) / np.sum(numcases)
avgw = np.sum(numevs * numcwooutl) / np.sum(numcwooutl)
poispred = np.asarray([np.sum(numcases) * pois(avg, x) for x in numevs])
poispredw = np.asarray([np.sum(numcases) * pois(avgw, x) for x in numevs])

print("Num: ", numevs)
print("Data: ", numcases)
print("Data w/o outliers: ", numcwooutl)
print("Poisson: ", poispred)
print("Poisson w/o o: ", poispred)
print("Average: ", avg)
print("Average w/o o: ", avgw)

plt.bar(numevs, numcases, width = 0.25, color='b')
plt.bar(numevs + 0.25, numcwooutl, width = 0.25, color='r')
plt.bar(numevs + 0.5, poispred, width = 0.25, color='g')
plt.bar(numevs + 0.75, poispredw, width = 0.25, color='k')


plt.xlabel(r"num events")
plt.ylabel(r"num intervals")

plt.savefig("ex8.2.eps", dpi=800)


#calculate chi2

def chi2(x, y):
    z = (x-y)**2 / y
    return np.sum(z)

print("Chi2 = ", chi2(numcases, poispred))
print("Chi2 w/o = ", chi2(numcwooutl, poispredw))
print("DoF = ", len(numcases) - 1)

print(f"Probability with outlier: {1 - scipy.stats.chi2.cdf(chi2(numcases, poispred), len(numcases) - 1) : .6e}")
print(f"1 - probability with outlier: {scipy.stats.chi2.cdf(chi2(numcases, poispred), len(numcases) - 1): .3E}")
print(f"Probability with outlier: {1 - scipy.stats.chi2.cdf(chi2(numcwooutl, poispredw), len(numcases) - 2) : .6e}")
print(f"1 - probability with outlier: {scipy.stats.chi2.cdf(chi2(numcwooutl, poispredw), len(numcases) - 2): .3E}")
