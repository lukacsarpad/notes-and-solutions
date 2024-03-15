#!/usr/bin/env python3

import numpy as np
from math import sqrt

marks = [[22, 63], [48, 39], [76, 61], [10, 30], [22, 51], [4, 44], [68, 74], [44, 78], [10, 55], [76, 58], [14, 41], [56, 69]]

print("Number of data points: ", len(marks))


def mean(x):
    mean = np.zeros(len(x[0]))
    n = len(x)
    for i in range(n):
        mean = mean + x[i]
    mean = mean / n
    return mean

print("Mean: ", mean(marks))

def stddev(x):
    stddev = np.zeros(len(x[0]))
    n = len(x)
    xbar = mean(x)
    for i in range(n):
        stddev = stddev + (x[i] - xbar) ** 2
    stddev = stddev / n
    stddev = np.sqrt(stddev)
    return stddev

print("Standard deviation: ", stddev(marks))

def cov(x):
    means = mean(x)
    n = len(x)
    m = len(x[0])
    covs = np.zeros([m, m])
    for nd in range(n):
        for i in range(m):
            for j in range(m):
                covs[i][j] = covs[i][j] + (x[nd][i] - means[i]) * (x[nd][j] - means[j])
    covs = covs/n
    return covs

print("Covariance matrix:\n", cov(marks))

def corr(x):
    covs = cov(x)
    n = len(covs)
    sigmas = [covs[i][i] for i in range(n)]
    sigmas = np.sqrt(sigmas)
    for i in range(n):
        for j in range(n):
            covs[i][j] = covs[i][j]/sigmas[i]/sigmas[j]
    return covs

print("Correlation matrix:\n", corr(marks))
