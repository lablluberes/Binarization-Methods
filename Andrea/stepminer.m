
%x1 - xn, n is time points (vector of increasing values)
%x is values (vector)
%X is resulting vector

function [u, threshold] = stepminer(x, alpha)

n = length(x);
%initialize
%setting is twostep function true or false
[SSR, SSE, t1, t2] = SSCalculate(n,x);
P1 = pval(SSE(1), SSR(1), n, 3);
P2 = pval(SSE(2), SSR(2), n, 4);
P12 = F12(SSE(1),SSE(2), 3, 4, n);

%determine which to use better

if P1<alpha % && (P12 > alpha)
    %type = "oneStep";
    u = x > t1;
    threshold = t1
elseif P2<alpha
    %type = "twoStep";
    u = x > t2;
    threshold = t2
else
    %type = "other";
    u = zeros(1,n);
    threshold = NaN
    
end


end


function [SSR, SSE, threshold1, threshold2] = SSCalculate(n, x)  %n here is the length not the vector

    %Calculate SSTOT

    xmean = mean(x);

    SSTOT = sum((x(1:n)-xmean).^2);

    %find means from left and right as it goes
    %and find smallest SSE for optimal step
    %considering all possibilities

    %rightMeanMin = rightMean
    %leftMeanMin = leftMean
    %breakIndex = -1

    SSEmin1 = SSTOT;
    SSEmin2 = SSTOT;

    %SSR = SSTOT - SSE
    %SSE = SSTOT - SSR
    %find min SSE

    for i = 1:n-1
        %add up mean from left side
        leftMean = mean(x(1:i));
        %add up mean from right side
        % this is so it doesnt divide by 0
        rightMean = mean(x(i+1:n));
        if i < n-1
            %all values from breakpoint i until possible breakpoints
            for j = i+1:n-1

                rightMean2 = mean(x(i+1:j));
                leftMean2 = (mean(x(1:i))*(i) + mean(x(j+1:n)) * (n-j))/(n-j+i);

                %first amount before breakpoint * leftMean
                %amount since breakpoint1 to breakpoint2
                %amount from breakpoint2 to end
                SSE2 = sum((x(1:i)-leftMean2).^2) + sum((x(i+1:j)-rightMean2).^2) + sum((x(j+1:n)-leftMean2).^2);

                if SSEmin2 > SSE2
                    SSEmin2 = SSE2;

                    %binarized vector
                    %values

                    threshold2 = mean([leftMean2 rightMean2]);

                    %breakpoints2(1) = leftMean2
                    %breakpoints2(2) = rightMean2
                    %breakpoints2(3) = i
                    %breakpoints2(4) = j

                end
                
            end


        %calculate SSE
        %calculate SSR first since xmean is fixed
        %and fixed values are constant until break point
    

        SSE1 = sum((x(1:i)-leftMean).^2) + sum((x(i+1:n)-rightMean).^2);
   
        %try every break point until min SSE
        if SSEmin1 > SSE1
            SSEmin1 = SSE1;

            %store means and breakpoints

            threshold1 = mean([leftMean rightMean]);

            %breakpoints1(1) = leftMean
            %breakpoints1(2) = rightMean
            %breakpoints1(3) = i
    
        end
    end

    %last return values
    SSE1 = SSEmin1;
    SSR1 = SSTOT - SSE1;

    SSE2 = SSEmin2;
    SSR2 = SSTOT - SSE2;

    SSR = [SSR1, SSR2];
    SSE = [SSE1, SSE2];

end


%F

function [P, m] = pval (SSE, SSR, n, m)



    MSR = SSR/(m-1);

    %MSE

    MSE = SSE/(n-m);

    %F

    F = MSR/MSE;

    %f distribution probability
    % i think this is correct
    %p value

    P = fpdf(F,m-1,n-m);


end

%F12

function F = F12(SSE1, SSE2, m1, m2, n)
    F = abs(((SSE1-SSE2)/(m2-m1))/(SSE2/(n-m2)));
    F = chi2cdf(F,1);
end



