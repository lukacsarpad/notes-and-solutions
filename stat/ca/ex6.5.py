#!/usr/bin/env python3

import numpy as np

times = np.asarray([1.0 * i for i in range(0, 9)])
counts = np.asarray([997, 520, 265, 127, 70, 35, 16, 7, 3])

# let the initial activity of the sample be A, then the activity at time t_i is exp(- log(2) * t / tau)
# yielding log(c_i) = log(A) - log(2) * t_i / tau
# let us solve for the inverse
# c_i = A * exp(-log(2) t / tau)
# log(c) = log(A) - log(2) t / tau
# t =  - tau /  log(2) * log(c) + tau / log(2) * log(A)

# weighting: there is a 1/sqrt(N) error on the counts
weights = counts ** 1.5
den = np.sum(weights)

x = times
y = np.log(counts)

xybar = np.sum(x * y * weights) / den
xbar = np.sum(x * weights) / den
ybar = np.sum(y * weights) / den
x2bar = np.sum(x * x * weights) / den

m = (xybar - xbar * ybar) / (x2bar - xbar**2)

c = ybar - m * xbar

tau = - np.log(2.) / m
act0 = np.exp(c)

# error estimate
n = len(times)
sigma2bar = n / den

sigma2m = sigma2bar / n / (x2bar - xbar**2)
errm = np.sqrt(sigma2m)
errtau = np.log(2) / m**2 * errm


print(f'Least squares fit: A = {act0 : .3f}, tau = {tau: .3f}.')
print(f'Error estimate for tau: {errtau: .3f}.')


    



