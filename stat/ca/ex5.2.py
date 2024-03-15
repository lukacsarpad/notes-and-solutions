#!/usr/bin/env python3

import numpy as np

nums = np.asarray([20.0, 19.7, 20.6, 18.5, 21.2, 20.8, 20.7])

mu = np.sum(nums)/len(nums)

print("Mean: ", mu)

devs = (nums - mu)**2

var = np.sum(devs) / (len(nums)-1)

print("Stddev: ", np.sqrt(var))

print("Err. estimate: ", 0.8/np.sqrt(7))

known_mean = 20.0
devs_known_mean =  (nums - known_mean)**2
var_known_mean  = np.sum(devs_known_mean) / len(nums)

print("Stddev with known mean: ", np.sqrt(var_known_mean))

print("Estimate of error of stddev w/o prior knowledge: ", np.sqrt(var/2/(len(nums)-1)))

print("Estimate on mean w. above: ", np.sqrt(var/7))




