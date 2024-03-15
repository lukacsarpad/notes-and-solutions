#!/usr/bin/env python3

import numpy as np
import math


resa=np.asarray([910, 950, 1050, 1060, 940, 1070, 1090, 930, 910, 1060])
resb=np.asarray([960, 1040, 980, 1010, 1000, 1020, 990, 990, 1010, 970])


# runtest
print("Performing run test.")
resa_wl = [[x, 'a'] for x in resa]
resb_wl = [[x, 'b'] for x in resb]

res_wl = resa_wl + resb_wl


def choose0(x):
    return x[0]

res_s = sorted(res_wl, key=choose0)

letts = [x[1] for x in res_s]

print("Letter sequence: ", ' '.join(letts))

runs = 1
for i in range(1, len(res_s)):
    if res_s[i][1] != res_s[i - 1][1]:
        runs = runs + 1
print("Number of runs: ", runs)

na = len(resa)
nb = len(resb)

expruns = 1 + 2 * na * nb / (na + nb)
varruns = 2 * na * nb*(2 * na * nb - na - nb) / (na + nb)**2 / (na + nb -1)

print("Expected number of runs: ", expruns)
print("Variance of runs: ", varruns)
print("Deviation: ", runs - expruns)
print("Divided by sigma: ", (runs - expruns) / np.sqrt(varruns))

print("Performing the Mann-Whitney test.")

ua = 0
ub = 0
for i in range(0, len(res_s)):
    ri = 0
    for j in range(i + 1, len(res_s)):
        if res_s[j][1] != res_s[i][1]:
            ri = ri + 1
    if res_s[i][1] == 'a':
        ua = ua + ri
    else:
        ub = ub + ri

print("Ua: ", ua, ", Ub: ", ub)

ra = 0
rb = 0 
for i in range(0, len(res_s)):
    if res_s[i][1] == 'a':
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



