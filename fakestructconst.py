# -*- coding: utf-8 -*-
"""
Created on Sat May  1 15:48:18 2021

@author: soric
"""
import numpy as np

def Inverse(i, n):
    return (i + int(round(n/2)))%n

def GenerateFakeStructconst(n, seed=None, r=2):
    if (type(n) != int) or (type(r) != int):
        raise Exception('Error: n and r must be integers')
    if (n<=0) or (r<=0):
        raise Exception('Error: n and r must be greater than zero')
    if (r%2 != 0):
        raise Exception('Error: r must be divisible by two')
    if (round((n - r)/2) < r):
        raise Exception('Error: (n-r)/2 must be greater than r')
    
    if type(seed) == int:
        if seed >= 0:
            rng = np.random.default_rng(seed)
        else:
            rng = np.random.default_rng()
    else:
        rng = np.random.default_rng()

    #There are (n-r)/2 positive roots and (n-r)/2 negative roots
    #The first r of the positive roots are the simple roots
    structconsts = np.zeros([n, n, n])
    for i in range(int(round(n/2))):
        for j in range(n):
            if max(structconsts[i, j] != 0) or (min(structconsts[i, j]) != 0):
                continue
            if j <= i:
                continue
            if i == j:
                continue
            nonzeroindex = rng.integers(low=0, high=n)
            nonzerovalue = rng.integers(low=0, high=4)*(-1)**(rng.integers(low=0, high=2))
            structconsts[i,j,nonzeroindex] = nonzerovalue
            structconsts[j,i, nonzeroindex] = -nonzerovalue
            Ini = Inverse(i, n)
            Inj = Inverse(j, n)
            Inz = Inverse(nonzeroindex, n)
            structconsts[Inj, Ini, Inz] = nonzerovalue
            structconsts[Ini, Inj, Inz] = -nonzerovalue

    return structconsts


