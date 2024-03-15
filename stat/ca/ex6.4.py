#!/usr/bin/env python3

import numpy as np

xv = np.asarray([0.1 * i for i in range(2, 8)])

yv = np.asarray([0.599, 0.896, 1.189, 1.479, 1.756, 2.044])


def chi2(x, y, a, b):
    yth = a * x + b * np.sin(x)
    return np.sum( (y - yth) ** 2)


# least squares fit
def fitlinsin(x, y):
    n = len(x)
    x2bar = np.sum(x * x) / n
    xsinxbar = np.sum(x * np.sin(x)) / n
    sin2xbar = np.sum(np.sin(x) ** 2) / n
#
    xybar = np.sum(x * y) / n
    ysinxbar = np.sum(np.sin(x) * y) / n
#
    det = x2bar * sin2xbar - xsinxbar ** 2
    a = (1 / det) * (sin2xbar * xybar - xsinxbar * ysinxbar)
    b = (1 / det) * (-xsinxbar * xybar + x2bar * ysinxbar)
    return a, b


a, b = fitlinsin(xv, yv)

print(f'Fitted value of a: {a : .3f}, and b: {b: .3f}')

    



