#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# TeX labels
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


x = [90, 90, 79, 84, 78, 91, 88, 90, 85, 80, \
     88, 75, 73, 79, 78, 79, 67, 83, 68, 60, \
     73, 79, 69, 74, 76, 68, 72, 72, 75, 60, \
     61, 66, 66, 54, 71, 67, 75, 49, 51, 57, \
     62, 64, 68, 58, 56, 79, 63, 68, 64, 51, \
     58, 53, 65, 57, 59, 65, 48, 54, 55, 40, \
     49, 42, 36, 46, 40, 37, 53, 48, 44, 43, \
     35, 39, 30, 41, 41, 22, 28, 36, 39, 51]

bins = [10 * i for i in range(11)]
hist, binEdges = np.histogram(x, bins)

print(hist)
print(binEdges)

plt.bar(binEdges[:-1], hist, width=np.diff(binEdges))

plt.xlabel(r"bins")
plt.ylabel(r"frequency")

plt.savefig("ex2.6.eps", dpi=800)

def mean(x):
    mean = 0
    n = len(x)
    for i in range(n):
        mean = mean + x[i]
    mean = mean / n
    return mean

def median(x):
    y = x
    y.sort()
    n = len(x)
    if((n % 2) == 0):
        return((y[n//2 - 1] + y[n//2])/2)
    else:
        return y[(n - 1) // 2]

print("Mean of data: ", mean(x))
print("Median of data: ", median(x))

def mode(x):
    y = x
    y.sort()
    minx = y[0]
    maxx = y[-1]
    numbins = maxx - minx + 1
    h, e = np.histogram(x, numbins)
    i = np.argmax(h)
    return y[0] + i

print("Mode of data: ", mode(x))

def mean_from_hist(h, e):
    n = len(h)
    nn = np.sum(h)
    avg = 0
    for i in range(n):
        y = (e[i] + e[i+1])/2
        print(y, h[i])
        avg = avg + y * h[i]
    avg = avg / nn
    return avg

# y = x
# y.sort()
# minx = y[0]
# maxx = y[-1]
# bins = range(minx, maxx+2)
# hist, binEdges = np.histogram(x, bins)

print("Average from histogram: ", mean_from_hist(hist, binEdges))

def med_from_hist(h, e):
    nn = np.sum(h)
    k = 0
    for i in range(len(h)):
        k = k + h[i]
        if k >= nn / 2:
            return (e[i] + e[i+1])/2

print("Median from histogram: ", med_from_hist(hist, binEdges))

imax = np.argmax(hist)
print("Mode from histogram: ", (binEdges[imax]+binEdges[imax+1])/2)

def stddev(x):
    stddev = 0
    n = len(x)
    xbar = mean(x)
    for i in range(n):
        stddev = stddev + (x[i] - xbar) ** 2
    stddev = stddev / n
    stddev = np.sqrt(stddev)
    return stddev

print("Standard deviation of the data: ", stddev(x))

# calculate fwhm
def xup(h, e):
    for i in range(imax, len(hist)):
        print("hi", h[i], h[imax]/2)
        if(h[i] < h[imax]/2):
            return(e[i])
    return e[-1]

def xdn(h, e):
    for i in range(imax, 0, -1):
        if(h[i] < h[imax]/2):
            return(e[i+1])
    return e[0]

print("FWHM: ", xup(hist, binEdges) - xdn(hist, binEdges))
print("Gassian FWHM: ", 2.35 * stddev(x))
