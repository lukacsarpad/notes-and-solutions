#!/usr/bin/env python3

import numpy as np

times = np.asarray([0.16, 0.40, 0.58, 0.72, 0.97])

distances = np.asarray([0.2, 1.0, 2.0, 3.0, 5.0])

sigma = 0.01

sqdists = np.sqrt(distances)

# dist = 1/2 g t^2
# t = sqrt(2 dist / g)

def chi2(sd, t, g, sigma):
    tth  = np.sqrt(2 / g) * sd
    sigma2 = sigma**2
    squares = (t - tth) **2 / sigma2
    return np.sum(squares)


def inv2sqgfit(t, sd):
    return np.sum(t * sd) / np.sum(sd ** 2)

g = 2 / inv2sqgfit(times, sqdists)**2

print("Least squares fit g: ", g)

print("chi^2: ", chi2(sqdists, times, g, sigma))

print("sqrt(2*chi2): ", np.sqrt(2*chi2(sqdists, times, g, sigma)))

print("DoF: ", np.sqrt(2*len(times) - 2))

def errestim(sd):
    return np.sqrt(sigma**2/np.sum(sd**2))

print("Error estimate on 1/sqrt(g): ", errestim(sqdists))

print("Error estimate on g: ", 4 * errestim(distances) / inv2sqgfit(times, sqdists) ** 3)

def fitnzi(t, sd):
    n = len(t)
    xybar = np.sum(t * sd) / n
    ybar = np.sum(t) / n
    xbar = np.sum(sd) / n
    x2bar = np.sum(sd ** 2) / n
    m = (xybar - xbar * ybar) / (x2bar - xbar ** 2)
    c = ybar - m * xbar
    return m, c

invsqg, t0 = fitnzi(times, sqdists)

print("Fitting with nonzero intercept: sqrt(2/g) = ", invsqg, ", t0 = ", t0)
print("Enhanced fit of g = ", 2 / invsqg**2)

def errestimnzig(t, sd, sigma):
    sigma2 = sigma ** 2
    x2bar = np.sum(sd**2) / len(sd)
    xbar = np.sum(sd) / len(sd)
    sigma2m = sigma2 / len(sd) / (x2bar - xbar**2)
    sigma2c = sigma2 * x2bar / len(sd) / (x2bar - xbar**2)
    rho = - xbar / np.sqrt(x2bar)
    return sigma2m, sigma2c, rho


vsig, vt0, corr = errestimnzig(times, sqdists, sigma)

print("Error estimate on g: ", 4 * np.sqrt(vsig) / invsqg ** 3)

def chi2nzi(sd, t, g, t0, sigma):
    tth = t0 + np.sqrt(2 / g) * sd
    sigma2 = sigma ** 2
    squares = (t - tth) ** 2 / sigma2
    return np.sum(squares)


print("chi^2 = ", chi2nzi(sqdists, times, 2 / invsqg**2, t0, sigma))
print("sqrt(2*chi2): ", np.sqrt(2*chi2nzi(sqdists, times, 2 / invsqg **2, t0, sigma)))

