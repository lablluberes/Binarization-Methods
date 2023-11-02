# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:47:54 2023

@author: kashi

binarize data based on biggest jump
shmulevich algorithm
"""

import numpy as np


def binarize(x):
    
    n = len(x)
    s = np.sort(x)
    d = np.empty(n)
    
    for i in range(n-2):
        d[i] = s[i+1] - s[i]
    
    t = (s[n-1] - s[0])/(n-1)
    
    mn = s[n-1]
    index = 0
    
    for i in range(n-1):
        if d[i] > t and d[i] < mn:
            mn = d[i]
            index = i
            
    z = s[index + 1]
   
    
    return z


def readMat(matrix):
    data = np.loadtxt(matrix,delimiter=",")
    blank = np.empty([len(data), 1])
    
    
    for i in range(len(data)):
        t = binarize(data[i])
        blank[i] = t
        print(t)
        
    #change name of text here
    np.savetxt("thresholds.txt", blank)
    return blank


#change name of doc here
readMat("HIVIn(Matlab).csv")
