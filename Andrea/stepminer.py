#This stepminer works as the paper describes, it verifies if the F value given
#for each optimal step function is greater than the critical F for the degrees
#of freedom for the F distribution based on the number of elements in the 
#vector. If F is below the critical value for both then the final vector will
#just result in 0 as neither onestep nor twostep are enough to reject the
#null hypothesis.


import numpy as np
import scipy.stats

#FUNCTIONS


def stepminer(x, alpha):
    n = len(x)
    
    
    SSE, SSR, t1, t2 = SScalculate(n,x)
    
    
    t = None
    
    #print(SSE)
    #print(SSR)
    
    m = 3
    
    #onestep
    MSR = SSR[0]/(m-1)
    MSE = SSE[0]/(n-m)
    F = MSR/MSE
    
    if F > scipy.stats.f.ppf(q=1-alpha, dfn=m-1, dfd=n-m):
        t = t1
        #print("step 1")
    else:      
         m = 4
         
         #twostep
         MSR = SSR[1]/(m-1)
         MSE = SSE[1]/(n-m)
         F = MSR/MSE
         if F > scipy.stats.f.ppf(q=1-alpha, dfn=m-1, dfd=n-m):
             t = t2
             #("step 2")
         #else:
            #print("no reg")
             

    #initialize at 0
             
    u = np.zeros(n)
    
    #change above threshold to 1's
    if t:
        for i in range(n):
            if x[i] >= t:
                u[i] = 1
          
    #return vector and threshold
    return u, t
    
def SScalculate(n, x):
    
    xmean = np.mean(x)
    SSTOT = getSSTOT(x, n, xmean)
    
    SSEmin1 = SSTOT
    SSEmin2 = SSTOT
    
    
    for i in range(n - 1):
        leftMean = np.mean(x[0:i+1])
    
        rightMean = np.mean(x[i+1:n])
        
        if i < n-2:
            for j in range(i+1, n-1):
                
                #verify later
                rightMean2 = np.mean(x[i+1:j+1])
                leftMean2 = (leftMean*(i+1) + np.mean(x[j+1:n])*(n-j-1))/(n-j+i)
                
                
                SSE2 = 0
                
                for k in range(n):
                    if k < i+1 or k > j:
                        SSE2 = SSE2 + (x[k]-leftMean2)**2
                    else:
                        SSE2 = SSE2 + (x[k] - rightMean2)**2
                        
                if SSEmin2 > SSE2:
                    SSEmin2 = SSE2
                    #print("2:",SSEmin2)
                    
                    t2 = (leftMean2 + rightMean2)/2
                    
        SSE1 = 0
        for j in range(n):
            if j < i+1:
                SSE1 = SSE1 + (x[j] - leftMean)**2
            else:
                SSE1 = SSE1 + (x[j] - rightMean)**2
                
        if SSEmin1 > SSE1:
            SSEmin1 = SSE1
            #print("1:",SSEmin1)
            
            t1 = (leftMean + rightMean)/2
            
            
    
    SSE = np.array([SSEmin1, SSEmin2])
    SSR = np.array([SSTOT - SSEmin1, SSTOT - SSEmin2])
    
    return SSE, SSR, t1, t2





    
def getSSTOT(x, n, xmean):
    m = 0
    for i in range(n):
        m = m + (x[i] - xmean)**2
    return m


#TESTING

HIVdat = np.loadtxt("HIVIn(Matlab).csv",delimiter=",")
blank = np.empty([len(HIVdat), len(HIVdat[0])])

for i in range(len(HIVdat)):
    arr, t = stepminer(HIVdat[i], 0.05)
    blank[i] = arr

np.savetxt("stepminerpy.txt", blank, delimiter=",", fmt='%d')