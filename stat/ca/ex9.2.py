#!/usr/bin/env python3

import numpy as np
import math

chipsa=np.asarray([41.2, 40.8, 38.8, 42.0, 41.6, 38.4, 39.8, 41.1, 38.1])
chipsb=np.asarray([38.9, 38.5, 40.5, 37.3, 36.9, 40.4, 40.6, 37.6])


# runtest
print("Performing run test.")
chipsa_wl = [[x, 'a'] for x in chipsa]
chipsb_wl = [[x, 'b'] for x in chipsb]

chips_wl = chipsa_wl + chipsb_wl


def choose0(x):
    return x[0]

chips_s = sorted(chips_wl, key=choose0)

letts = [x[1] for x in chips_s]

print("Letter sequence: ", ' '.join(letts))

runs = 1
for i in range(1, len(chips_s)):
    if chips_s[i][1] != chips_s[i - 1][1]:
        runs = runs + 1
print("Number of runs: ", runs)

na = len(chipsa)
nb = len(chipsb)

expruns = 1 + 2 * na * nb / (na + nb)
varruns = 2 * na * nb*(2 * na * nb - na - nb) / (na + nb)**2 / (na + nb -1)

print("Expected number of runs: ", expruns)
print("Variance of runs: ", varruns)
print("Deviation: ", runs - expruns)
print("Divided by sigma: ", (runs - expruns) / np.sqrt(varruns))

print("Performing the Mann-Whitney test.")

ua = 0
ub = 0
for i in range(0, len(chips_s)):
    ri = 0
    for j in range(i + 1, len(chips_s)):
        if chips_s[j][1] != chips_s[i][1]:
            ri = ri + 1
    if chips_s[i][1] == 'a':
        ua = ua + ri
    else:
        ub = ub + ri

print("Ua: ", ua, ", Ub: ", ub)

ra = 0
rb = 0 
for i in range(0, len(chips_s)):
    if chips_s[i][1] == 'a':
        ra = ra + (i + 1)
    else:
        rb = rb + (i + 1)

print("Ra: ", ra, ", Rb: ", rb)
print("Ua: ", round(na * nb + na * (na + 1)/2 - ra), ", Ub: ", round(na * nb + nb * (nb + 1)/2 - rb))

print("na * nb: ", na * nb, ", ua + ub: ", ua + ub)

print("Expected Ua: ", na * nb / 2)
print("Variance: ", 1/12 * na * nb * (na + nb))
print("Deviation: ", ua - na * nb / 2)
print("Divided by sigma: ", (ua - na * nb / 2) / np.sqrt(1/12 * na * nb * (na + nb)))


