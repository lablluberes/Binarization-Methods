% Stepminer algoritm 

% Not finished

% To run use Stepminer(genes = vector) the function returns means (not
% finished) 
function step = Stepminer(genes)

    %genes = sort(genes);

    % length of genes vector
    n = length(genes);

    % calculate SSTOT
    SSTOT = sum((genes - mean(genes)).^2);

    minSSE1 = Inf;
    minSSE2 = Inf;

    % for each element in the vector calculate the fitted value
    % based on mean from point 1 to i on left side and mean of i+1 to n on
    % the right side
    for i = 1:n
        % left side mean for fitted value
        meansOneLeft = mean(genes(1:i));

        % right side mean for fitted value
        if (i+1 <= n)
            meansOneRight = mean(genes(i+1:n));
        end

        % find twostep sse and mean values
        for j = i+1:n
            % mean on left side
            meansTwoLeft = mean(genes(1:i));

            % mean of the middle
            %if(j+1 <= n)
            meansTwoMiddle = mean(genes(i+1:j));
            %end

            % mean of right side
            meansTwoRight = mean(genes(j+1:n));

            % calculate sse2
            SSE2 = sum(((genes(1:i) - meansTwoLeft).^2)) + sum(((genes(i+1:j) - meansTwoMiddle).^2)) + sum(((genes(j+1:n) - meansTwoRight).^2));
          
            % get the minimum sse2
            if (minSSE2 > SSE2)
                minSSE2 = SSE2;
            end

        end


        % calculate the sse of those means (sum it in order to look for the
        % minimum value)
        SSE1 = sum(((genes(1:i) - meansOneLeft).^2)) + sum(((genes(i+1:n) - meansOneRight).^2));

        % get the minium sse1
        if (minSSE1 > SSE1)
            minSSE1 = SSE1;
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
    F12 = ((minSSE1 - minSSE2) / (4 - 3)) / (minSSE2 / (n - 4));

    % p-values
    P1 = fpdf(F1,3-1,n-3);
    P2 = fpdf(F1,4-1,n-4);

    % decide if it is one, two step or other
    if (P1 < 0.05 && (F12 >= 0.05))
        step = "OneStep";

    elseif (P2 < 0.05)
        step = "TwoStep";
    
    else
        step = "Other";
    end

end