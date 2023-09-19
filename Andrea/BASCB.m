
%returns binary data



function u = BASCB(x) %x is vector of time series measurements, sorted in
                        %ascending order
    %initial step function is already given with vector of ascending order
    

    %sort vector
    
    sorted = sort(x)

    %step 1

    %length of array
    n = length(x)
    
    %default sigma
    sigma=[0.01:0.001:0.02]
    %some functions

    %I = besselj(n,x)

    %delta
    %delta(i) = x(i+1) - x(i) 



    %delta sigma
    %vector of length n
    %to store deltaSigma
    %n-1 so it doesnt try to access index n+1
    deltaSigma = [1:n-1]
    %vector of vectors, 2d vector
    deltaSigmaVec = zeros(length(sigma),n-1)
    maxSigmaVec = zeros(length(sigma),n-1)
    %for maxima
    maxSigma = [1:n-1]

    %for every sigma parameter
    for sig = 1:length(sigma)
        %try smoothing every step
        for i = 1:n-1
            sigsum = 0
    
            for j = -10:10
            %simulate -inf to inf???
    
                sigsum = sigsum + ((sorted(i+1) - sorted(i))*T(i-j,sig))
    
            end
            
            deltaSigma(i) = sigsum
        end
       %assign vector to location on vector of vectors

       for i = 1:n-1
            deltaSigmaVec(sig,i) = deltaSigma(i)
       end

       %local maxima

       for j = 1:n-1
    
          res = 0
          if j == 1
              res = deltaSigmaVec(sig,j+1) < deltaSigmaVec(sig,j)
          elseif j == n-1
              res = deltaSigmaVec(sig,j-1) < deltaSigmaVec(sig,j)
          else
              res = deltaSigmaVec(sig,j-1) < deltaSigmaVec(sig,j) & deltaSigmaVec(sig, j+1) < deltaSigmaVec(sig,j)
          end
          %add 1's and 0's to maxSigma array
          
          maxSigma(j) = res

       end


       %assign values
       for i = 1:n-1
            maxSigmaVec(sig,i) = maxSigma(i)
       end

    end


    %traverse maxSigmaVec aka Msigma. find funtion with single remaining
    %discontinuity ????
    %i dont know what that means

    
    %delete doubles

    [MaxSig, MaxSigIndex] = unique(maxSigmaVec,'rows')

    %step 2
    
    v = zeros(1,length(transpose(MaxSigIndex)))
    counter = 1
    
    for k = transpose(MaxSigIndex)
    smoothed = 1:n-1
    qVec = zeros(1,n-1)
    

       %smoothed function
       for i = 1:n
    
            smoothsum = 0
            for j = 1:i-1
                smoothsum = smoothsum + deltaSigmaVec(k,j)
            end
            %u1 isfirst value of
            %original function
            u1 = x(1)
            smoothsum = smoothsum + u1 
            smoothed(i) = smoothsum
        end



            %discontinuity h

            j = sum(maxSigmaVec(k,:))
            breakpoints = find(maxSigmaVec(k,:)==1) 

            for i = 1:j

                if i == 1 & j > 1
                    h = Y_ab(smoothed,breakpoints(i)+1,breakpoints(i+1)) - Y_ab(smoothed,1,breakpoints(i))
                elseif i == j & j>1
                    h = Y_ab(smoothed,breakpoints(i)+1,n) - Y_ab(smoothed,breakpoints(i-1)+1,breakpoints(i))
                elseif i == j & j == 1
                    h = Y_ab(smoothed,breakpoints(i)+1,n) - Y_ab(smoothed,1,breakpoints(i))
                else
                    h = Y_ab(smoothed,breakpoints(i)+1,breakpoints(i+1)) - Y_ab(smoothed,breakpoints(i-1)+1,breakpoints(i))
                end


                z = (sorted(breakpoints(i)) + sorted(breakpoints(i)+1))/2

                e = sum((sorted(i:n)-z).^2)

                q = h/e

                qVec(i) = q
            end

            [val, index] = max(qVec)
            v(counter) = index

            counter = counter + 1
   end
    

    
    %find v vector

    %[y, z] = max(qSigVec)
    %breakpoints = find(maxSigmaVec(z,:)==1) 
    %max of each col
    %y is the max val
    %z is the index
    
    %array of indexes basically
    %every spot where there could be a breakpoint
    %v = breakpoints

    %end

    
    %step 3
    %estimate location and variation of strongest discontinuities

    vmed = median(v)
    vmed = round(vmed)
    
    t = (sorted(vmed+1)+sorted(vmed))/2

    %binarized vector calculation

    %init binarized vector
    u = [1:n]

    for i = 1:n
           
        if x(i) <= t
            u(i) = 0
        else
            u(i) = 1
        end      
    end

end

function z = T(n,sig)
    %sigma is the smoothing parameter?
    z = exp(-2*sig)*besselj(n,2*sig)
end

function val = Y_ab(x,a,b)

    val = sum(x(a:b))/(b-a+1)

end