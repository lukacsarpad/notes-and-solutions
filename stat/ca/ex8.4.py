#!/usr/bin/env python3

import numpy as np
import math

temps = np.asarray([10.2, 10.4, 9.8, 10.5, 9.9, 9.8, 10.3, 10.1, 10.3, 9.9])
prec = 0.2

avg  = np.sum(temps) / len(temps)

devs = (temps - avg) ** 2 / prec**2

chi2 = np.sum(devs)

print(f"Average: {avg: .4f} +- {0.2 / np.sqrt(len(temps)) : .4f}")
print(f"chi2: {chi2: .4f}")
print(f"sqrt 2 chi2: {np.sqrt(2*chi2) : .4f}")
print(f"DoF: {len(temps) - 1}")

truet = 10.1
chi2a = np.sum((temps - truet) **2 / prec**2)

print(f"chi2 assuming true value of 10.1 K: {chi2a: .4f}")
print(f"DoF in this case: {len(temps)}")
print(f"sqrt 2 chi2: {np.sqrt(2*chi2a) : .4f}")
