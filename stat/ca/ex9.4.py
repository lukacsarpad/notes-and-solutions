#!/usr/bin/env python3

import numpy as np
import math
import sys


hb = np.asarray([80, 97, 110, 94, 120, 77, 84, 80, 87, 93, 120, 101, 82, 94, 100])
wf = np.asarray([85, 96, 101, 93, 122, 78, 84, 110, 79, 92, 110, 119, 91, 97, 111])

nn = len(hb)

if len(wf) != nn:
    print("Data input error.")
    sys.exit(-1)

# Pearson correlation
avgh = np.sum(hb) / nn
avgw = np.sum(wf) / nn

varh = np.sum((hb - avgh) ** 2) / nn
varw = np.sum((wf - avgw) ** 2) / nn

corrwh = np.dot((hb - avgh), (wf - avgw)) / nn

print(f"Number of couples: {nn}.")
print(f"Average of husbands: {avgh: .3f}, wifes: {avgw: .3f}.")
print(f"Variance of h: {varh: .3f}, w: {varw: .3f}.")
print(f"Stddev of h: {np.sqrt(varh): .3f}, w: {np.sqrt(varw): .3f}.")
print(f"The Pearson correlation coefficient is {corrwh / np.sqrt(varh * varw) : .3f}.")

# Spearman's correlation coefficient
hbn = [[i, hb[i]] for i in range(0, nn)]
wfn = [[i, wf[i]] for i in range(0, nn)]

def choose1(x):
    return x[1]


hbs = sorted(hbn, key=choose1)
wfs = sorted(wfn, key=choose1)

hbr = np.zeros(nn)
wfr = np.zeros(nn)

for i in range(0, nn):
    hbr[hbs[i][0]] = i
    wfr[wfs[i][0]] = i

diffsr = hbr - wfr
spearman = 1 - 6 * np.sum(diffsr ** 2) / (nn **3 - nn)
print(f"Spearman's correlation coefficient: {spearman: .3f}.")

diffs = hb - wf
diffsnums = [[i, diffs[i]] for i in range(0, nn)]
diffsnums = [x for x in diffsnums if x[1] != 0]
nnz = len(diffsnums)

def chooseabs1(x):
    return np.abs(x[1])

diffsnumss = sorted(diffsnums, key = chooseabs1)

diffsnumssr = np.zeros(len(diffsnumss))
i=0
while i < nnz:
    jj = [j for j in range(i, nnz) if np.abs(diffsnumss[j][1]) == np.abs(diffsnumss[i][1]) ]
    diffsnumssr[jj] = np.sum(np.asarray(jj) + 1)/len(jj)
    i = max(jj) + 1

ranksp = [diffsnumssr[i] for i in range(0, nnz) if diffsnumss[i][1] > 0]
ranksn = [diffsnumssr[i] for i in range(0, nnz) if diffsnumss[i][1] < 0]

srp = np.sum(ranksp)
srn = np.sum(ranksn)

print("Number of nonzero differences: ", nnz, ".")
print("Sum of ranks where difference positive: ", srp, ", where difference negative: ", srn, ", sum: ", srp + srn, ", n(n+1)/2: ", nnz * (nnz + 1) / 2, ".")
print("Smaller of the rank sums: ", min(srp, srn), ".")

