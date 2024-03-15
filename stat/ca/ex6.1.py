#!/usr/bin/env python3

import numpy as np

times = np.asarray([1.0*i for i in range(1, 7)])

distances = np.asarray([11., 19., 33., 40., 49., 61.])

sigma = 2.0

# Aim: least squares fit, with known intercept
# The chi2 to be minimised is

def chi2(t, u, v, sigma):
    uth = times * v
    sigma2 = sigma**2
    squares = (u - uth) **2 / sigma2
    return np.sum(squares)

# Normal equations: derivative of sigma w.r.t. v
# sum_i - 2 * (u_i - v * t_i) * t_i / sigma2
# solution is obtained by performing the summation as
# 2 * v * sum_i t_i**2 / sigma2 = 2 * sum_i u_i * t_i / sigma2
# divide by sigma2, as that is the same everywhere
# 2 * v * sum_i t_i**2 = 2 * sum_i  u_i t_i,
# v = (sum_i u_i t_i) / sum_i t_i**2

def vfit(t, u):
    return np.sum(u * t) / np.sum(t ** 2)

velocity = vfit(times, distances)

print("Least squares fit velocity: ", velocity)

print("chi^2: ", chi2(times, distances, velocity, sigma))

print("sqrt(2*chi2): ", np.sqrt(2*chi2(times, distances, velocity, sigma)))

print("DoF: ", np.sqrt(2*len(times) - 2))

def errestim(u):
    return np.sqrt(sigma**2/np.sum(u**2))

print("Error estimate on the velocity: ", errestim(distances))





