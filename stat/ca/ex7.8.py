#!/usr/bin/env python3

import numpy as np
from math import sqrt

pcs = np.asarray([6, 7, 9, 8])

print("Number of data points: ", len(pcs))


def mean(x):
    mean = 0
    n = len(x)
    for i in range(n):
        mean = mean + x[i]
    mean = mean / n
    return mean

print("Mean: ", mean(pcs))

def stddev(x):
    stddev = 0
    n = len(x)
    xbar = mean(x)
    for i in range(n):
        stddev = stddev + (x[i] - xbar) ** 2
    stddev = stddev / (n-1)
    stddev = sqrt(stddev)
    return stddev


print("Standard deviation: ", stddev(pcs))
err = stddev(pcs)/sqrt(len(pcs))
print("Error estimate for the mean: ", err)

mu=5.6
t = (mean(pcs) - mu) / err
print("t-value: ", t)

