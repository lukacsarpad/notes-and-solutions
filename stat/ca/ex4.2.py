#!/usr/bin/env python3

from math import sqrt

meas = [299798000, 299789000, 299797000]
errs = [5000, 4000, 8000]

def combine_meas(m, e):
    n = len(m)
    avg = 0.0
    denom = 0.0
    for i in range(n):
        avg = avg + m[i]/e[i]**2
        denom = denom + 1/e[i]**2
    avg = avg / denom
    return avg, 1/denom

cm, ce = combine_meas(meas, errs)


print("Combined measurement: ", cm, ", variance: ", ce, ", error: ", sqrt(ce))
