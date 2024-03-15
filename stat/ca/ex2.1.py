#!/usr/bin/env python3

from math import sqrt

ages=[19.0, 18.7, 19.3, 19.2, 18.9, 19.0, 20.2, 19.9, 18.6, 19.4, 19.3, 18.8, 19.3, 19.2, 18.7, 18.5, 18.6, 19.7, 19.9, 20.0, 19.5, 19.4, 19.6, 20.0, 18.9]

print("Number of data points: ", len(ages))


def mean(x):
    mean = 0
    n = len(x)
    for i in range(n):
        mean = mean + x[i]
    mean = mean / n
    return mean

print("Mean: ", mean(ages))

def stddev(x):
    stddev = 0
    n = len(x)
    xbar = mean(x)
    for i in range(n):
        stddev = stddev + (x[i] - xbar) ** 2
    stddev = stddev / n
    stddev = sqrt(stddev)
    return stddev


print("Standard deviation: ", stddev(ages))

def skew(x):
    skew = 0
    n = len(x)
    xbar = mean(x)
    sigma = stddev(x)
    for i in range(n):
        skew = skew + (x[i] - xbar)**3
    skew = skew / n
    skew = skew / sigma**3
    return skew

print("Skew: ", skew(ages))
