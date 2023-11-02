# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 10:29:51 2023

@author: kashi

This will return the binarization threshold + the minimum threshold so far
and the maximum threshold so far

Takes a matrix of data as input
"""
import numpy as np

def getSSTOT(x, n, xmean):
    m = 0
    for i in range(n):
        m = m + (x[i] - xmean)**2
    return m


def onestep(x):
    
    n = len(x)
    #step = 0
    xmean = np.mean(x)
    SSTOT = getSSTOT(x, n, xmean)
    
    SSEmin = SSTOT
    
    for i in range(n-1):
        leftMean = np.mean(x[0:i+1])
    
        rightMean = np.mean(x[i+1:n])
        
        SSE = 0
        
        for j in range(n):
            if j < i+1:
                SSE = SSE + (x[j] - leftMean)**2
            else:
                SSE = SSE + (x[j] - rightMean)**2
                    
        
        if SSEmin > SSE:
            SSEmin = SSE
            #print("1:",SSEmin1)
                
            t = (leftMean + rightMean)/2
    
    
    
    return t



def readMat(matrix):
    data = np.loadtxt(matrix,delimiter=",")
    blank = np.empty([len(data), 1])
    
    
    for i in range(len(data)):
        t = onestep(data[i])
        blank[i] = t
        #print(t)
        
    #change name of text here
    np.savetxt("thresholds.txt", blank)
    return blank



#replace this with name of matrix to read
readMat("HIVIn(Matlab).csv")