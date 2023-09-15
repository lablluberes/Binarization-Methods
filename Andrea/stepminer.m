
%x1 - xn, n is time points (vector of increasing values)
%x is values (vector)
%X is resulting vector

function [breakPoints, type] = stepminer(x, alpha)

n = length(x)
%initialize
%setting is twostep function true or false
[SSR, SSE, breakpoints1, breakpoints2] = SSCalculate(n,x)
[P1, m1] = pval(SSE(1), SSR(1), n)
[P2, m2] = pval(SSE(2), SSR(2), n)
P12 = F12(SSE(1),SSE(2), m1, m2,n)

%determine which to use better

if P1<alpha && P12 >= alpha
    type = "oneStep"
    breakPoints = breakpoints1
elseif P2<alpha
    type = "twoStep"
    breakPoints = breakpoints2
else
    type = "other"
    breakPoints = [mean(x), NaN, NaN, NaN]
    
end

type
breakPoints

end


function [SSR, SSE, breakpoints1, breakpoints2] = SSCalculate(n, x)  %n here is the length not the vector

    %Calculate SSTOT

    xmean = mean(x)

    SSTOT = 0
    for index = x
        SSTOT = SSTOT + (index - xmean)^2
    end


    %heavily based on 
    %https://github.com/suryattheja/StepMiner-BioInformatics

    leftMean = 0
    rightMean = xmean
    rightMean2 = mean(x(1:n-1)) %all but last value

    n1 = 1 %left start
    n2 = n %right start

    %index 1 is mean 1
    %index 2 is mean 2
    %index 3 is breakpoint 1 index
    %index 4 is breakpoint 2 index (unused so NaN)
    breakpoints1 = [0 0 0 NaN]
    breakpoints2 = [0 0 0 0]

    %find means from left and right as it goes
    %and find smallest SSE for optimal step
    %considering all possibilities

    %rightMeanMin = rightMean
    %leftMeanMin = leftMean
    %breakIndex = -1

    SSEmin1 = SSTOT
    SSE1 = SSTOT
    SSEmin2 = SSTOT
    SSE2 = SSTOT

    %SSR = SSTOT - SSE
    %SSE = SSTOT - SSR
    %find min SSE

    for i = 1:n
        %add up mean from left side
        leftMean = mean(x(1:i))
        %add up mean from right side
        % this is so it doesnt divide by 0
        if n1 < n-1
            rightMean = mean(x(i:n))
            %all values from breakpoint i until possible breakpoints
            for j = i:n-1

                rightMean2 = mean(x(i+1:j))

                SSR2 = (n1*((leftMean-xmean)^2)) + ((j-i)*((rightMean2 - xmean)^2)) + (leftMean - xmean)^2
                SSE2 = SSTOT - SSR2

                if SSEmin2 > SSE2
                SSEmin2 = SSE2

                %values

                breakpoints2(1) = leftMean
                breakpoints2(2) = rightMean2
                breakpoints2(3) = i
                breakpoints2(4) = j

                end
            end

        n2 = n2 - 1 
        elseif n1 < n
            rightMean = mean(x(i:n))
            n2 = n2 - 1 
        else
            n2 = 1
        end

        %calculate SSE
        %calculate SSR first since xmean is fixed
        %and fixed values are constant until break point
    

        SSR1 = (n1 * ((leftMean - xmean)^2)) + (n2 * ((rightMean - xmean)^2))
        SSE1 = SSTOT - SSR1
        
        n1 = n1 + 1
        %try every break point until min SSE
        if SSEmin1 > SSE1
            SSEmin1 = SSE1

            %store means and breakpoints

            breakpoints1(1) = leftMean
            breakpoints1(2) = rightMean
            breakpoints1(3) = i
    
        end
    end

    %last return values
    SSE1 = SSEmin1
    SSR1 = SSTOT - SSE1

    SSE2 = SSEmin2
    SSR2 = SSTOT - SSE2

    SSR = [SSR1, SSR2]
    SSE = [SSE1, SSE2]

end



%F

function [P, m] = pval (SSE, SSR, n)


    if n > 4
        m = 4
    elseif n < 3
        m = 1
    else
        m = 3
    end

    %MSR

    MSR = SSR/(m-1)

    %MSE

    MSE = SSE/(n-m)

    %F

    F = MSR/MSE

    %f distribution probability
    % i think this is correct
    %p value

    P = fpdf(F,m-1,n-m)


end

%F12

function P = F12(SSE1, SEE2, m1, m2, n)
    P = ((SSE1-SEE2)/(m2-m1))/(SEE2/(n-m2))
end



