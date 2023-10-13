#This stepminer is modified to give the optimal step function between onestep
#and twostep based on which F value is greater and threfore more significant.
#The resulting vector will always include at least one 1 or at least one 0, as
#each resulting vector will be the result of a step function. The downside to
#this one is, sometimes even the best F value is too low to be statistically
#significant, yet it'll still accomodate a function for it.



import numpy as np

#FUNCTIONS


def stepminer(x, alpha):
    n = len(x)
    #step = 0
    
    xmean = np.mean(x)
    SSTOT = getSSTOT(x, n, xmean)
    
    SSEmin1 = SSTOT
    SSEmin2 = SSTOT
    
    #init 
    u = np.zeros(n)
    maxF = float('-inf')
    
    
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

                    
                    #verify F best fit
                    SSR = SSTOT - SSEmin2
                    m = 4
                    MSR = SSR/(m-1)
                    MSE = SSEmin2/(n-m)
                    F = MSR/MSE
                    
                    if F > maxF:
                        maxF = F
                        if rightMean2 > leftMean2:                           
                            for k in range(n):
                                if k >= i+1 and k < j+1:
                                    u[k] = 1
                                else:
                                    u[k] = 0
                        else:
                            for k in range(n):
                                if k >= i+1 and k < j+1:
                                    u[k] = 0
                                else:
                                    u[k] = 1
                                
                        t = (leftMean2 + rightMean2)/2
                        #step = 2
                    
                    
                    
        SSE1 = 0
        for j in range(n):
            if j < i+1:
                SSE1 = SSE1 + (x[j] - leftMean)**2
            else:
                SSE1 = SSE1 + (x[j] - rightMean)**2
                
        if SSEmin1 > SSE1:
            SSEmin1 = SSE1
            
            #verify F best fit
            SSR = SSTOT - SSEmin1
            m = 3
            MSR = SSR/(m-1)
            MSE = SSEmin1/(n-m)
            F = MSR/MSE
                    
            if F > maxF:
                maxF = F
                if rightMean > leftMean:   
                    for k in range(n):
                        if k >= i+1:
                            u[k] = 1
                        else:
                            u[k] = 0
                else:
                    for k in range(n):
                        if k >= i+1:
                            u[k] = 0
                        else:
                            u[k] = 1
                                
                t = (leftMean + rightMean)/2
                #step = 1
                    
    
    #print("step=",step)
    #return vector and threshold
    return u, t

    
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

np.savetxt("stepminermodpy.txt", blank, delimiter=",", fmt='%d')