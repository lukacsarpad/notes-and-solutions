#!/usr/bin/env python3

import numpy as np

times = np.asarray([1.1, 2.2, 2.9, 4.1, 5.0, 5.8])

distances = np.asarray([10.0 * i for i in range(1, 7)])

sigma = 0.1

def chi2(u, t, v, sigma):
    tth  = u / v
    sigma2 = sigma**2
    squares = (t - tth) **2 / sigma2
    return np.sum(squares)


def vinvfit(t, u):
    return np.sum(u * t) / np.sum(u ** 2)

velocity = 1/vinvfit(times, distances)

print("Least squares fit velocity: ", velocity)

print("chi^2: ", chi2(distances, times, velocity, sigma))

print("sqrt(2*chi2): ", np.sqrt(2*chi2(distances, times, velocity, sigma)))

print("DoF: ", np.sqrt(2*len(times) - 2))

def errestim(d):
    return np.sqrt(sigma**2/np.sum(d**2))

print("Error estimate on 1/v: ", errestim(distances))

print("Error estimate on v: ", errestim(distances) / vinvfit(times, distances) ** 2)





