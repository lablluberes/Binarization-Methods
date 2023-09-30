% Stepminer algoritm 

% Not finished

% To run use Stepminer(genes = vector) the function returns means (not
% finished) 
function [binary, step, thr] = Stepminer(genes)

    %genes = sort(genes);

    % length of genes vector
    n = length(genes);

    % calculate SSTOT
    SSTOT = sum((genes - mean(genes)).^2);

    minSSE1 = Inf;
    minSSE2 = Inf;

    % for each element in the vector calculate the fitted value
    % based on mean from point 1 to i on left side and mean of i+1 to n on
    % the right side. Also calculate the middle means for two step from i+1
    % to j for the two step
    for i = 1:n-1
        % left side mean for fitted value
        meansOneLeft = mean(genes(1:i));

        % right side mean for fitted value
        %if (i+1 < n)
        meansOneRight = mean(genes(i+1:n));
        %end

        % find twostep sse and mean values
        for j = i+1:n-1
            % mean on left side and right side not middle
            meansTwoLeftRight = mean([genes(1:i) genes(j+1:n)]);

            % mean of the middle
            %if(j+1 <= n)
            meansTwoMiddle = mean(genes(i+1:j));
            %end

            % calculate sse2
            SSE2 = sum((genes(1:i) - meansTwoLeftRight).^2) + sum((genes(i+1:j) - meansTwoMiddle).^2) + sum((genes(j+1:n) - meansTwoLeftRight).^2);
          
            % get the minimum sse2
            if (minSSE2 > SSE2)
                minSSE2 = SSE2;

                % Save the means to a threshold vector to calculate it
                % later
                thresholdTwo(1:i) = meansTwoLeftRight;
                thresholdTwo(j+1:n) = meansTwoLeftRight;
                thresholdTwo(i+1:j) = meansTwoMiddle;
            end

        end


        % calculate the sse of those means (sum it in order to look for the
        % minimum value)
        SSE1 = sum((genes(1:i) - meansOneLeft).^2) + sum((genes(i+1:n) - meansOneRight).^2);

        % get the minium sse1
        if (minSSE1 > SSE1)
            minSSE1 = SSE1;

            % Save the means to a threshold vector to calculate it
            % later
            thresholdOne(1:i) = meansOneLeft;
            thresholdOne(i+1:n) = meansOneRight;
        end

    end


    % get the minimum sse values
    %[minSSE1, index1] = min(SSE1);
    %[minSSE2, index2] = min(SSE2);


    % calculate SSR for both one and two step
    SSR1 = SSTOT - minSSE1;
    SSR2 = SSTOT - minSSE2;

    % calculate MSR and MSE for one step
    MSR1 = SSR1 / (3 - 1);
    MSE1 = minSSE1 / (n - 3);
    
    % calculate MSR and MSE for two step
    MSR2 = SSR2 / (4 - 1);
    MSE2 = minSSE2 / (n - 4);

    % F statistic for one step
    F1 = MSR1 / MSE1;

    % F statistic for two step
    F2 = MSR2 / MSE2;

    % F statistic represents the relative goodness of fit of a one-step 
    % versus a two-step pattern.
    F12 = abs(((minSSE1 - minSSE2) / (4 - 3)) / (minSSE2 / (n - 4)));
    F12 = chi2cdf(F12, 1);

    % p-values
    P1 = fpdf(F1, 3-1, n-3);
    P2 = fpdf(F2, 4-1, n-4);

    % decide if it is one, two step or other
    if (P1 < 0.05 && (F12 >= 0.05))
        % assign step value
        step = "OneStep";

        % calculate threshold
        thr = mean(unique(thresholdOne));

        % look through vector if the thr is less than the mean then 0 or
        % else its 1
        for i = 1:n

            if (thresholdOne(i) <= thr)
                binary(i) = 0;
            else
                binary(i) = 1;
            end

        end

    elseif (P2 < 0.05)
        % assign step value
        step = "TwoStep";

        % calculate threshold
        thr = mean(unique(thresholdTwo));

        % look through vector if the thr is less than the mean then 0 or
        % else its 1
        for i = 1:n
            
            if (thresholdTwo(i) <= thr)
                binary(i) = 0;
            else
                binary(i) = 1;
            end

        end
    
    else
        % assign step value
        step = "Other";
        thr = 0;
        % assigned binary vector -1 because it is other
        for i = 1:n
            binary(i) = -1;
        end
    end

end

