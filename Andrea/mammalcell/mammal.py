# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 11:48:26 2023

@author: kashi
"""

import numpy as np


#indexes and what they represent
#0 = CycD
#1 = Rb
#2 = E2F
#3 = CycE
#4 = CycA
#5 = p27
#6 = Cdc20
#7 = Cdh1
#8 = UbcH10
#9 = CycB

def mammalcell(x):
    
    u = np.empty(10)
    
    u[0] = x[0]
    u[1] = not(x[0] or x[3] or x[4] or x[9]) or (x[5] and not(x[1] or x[9]))
    u[2] = not(x[1] or x[4] or x[9]) or (x[5] and not(x[1] or x[9]))
    u[3] = x[2] and not x[1]
    u[4] = (x[2] and not(x[1] or x[6]) and not (x[7] and x[8])) or (x[4] and not (x[1] or x[6]) and not(x[7] and x[8]))
    u[5] = not(x[0] or x[3] or x[4] or x[9]) or (x[5] and not(x[3] and x[4]) and not(x[9] or x[0])) 
    u[6] = x[9]
    u[7] = not(x[4] or x[9]) or x[6] or (x[5] and not x[9])
    u[8] = not(x[7]) or (x[7] and x[8] and (x[6] or x[4] or x[9]))
    u[9] = not(x[6] or x[7])
    
    
    return u




def readMat(matrix):
    data = np.loadtxt(matrix)
    blank = np.empty([len(data), 10])
    
    
    for i in range(len(data)):
        t = mammalcell(data[i])
        blank[i] = t
        #print(t)
        
    #change name of text here
    np.savetxt("mammalshmulevich.txt", blank, fmt='%d')
    return blank



readMat("shmulevich.txt")