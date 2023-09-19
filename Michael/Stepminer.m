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

    % for each element in the vector calculate the fitted value
    % based on mean from point 1 to i on left side and mean of i+1 to n on
    % the right side
    for i = 1:n
        % left side mean for fitted value
        means(i, 1) = mean(genes(1:i));

        % right side mean for fitted value
        if (i+1 <= n)
            means(i, 2) = mean(genes(i+1:n));
        end

        % find twostep sse and mean values
        for j = i:n
            % mean on left side
            means2(j, 1) = mean(genes(i:j));

            % mean of right side
            if(j+1 <= n)
                means2(j, 2) = mean(mean(genes(j+1:n)) + mean(genes(1:j)));
            end

            % calculate sse2
            SSE2(i, 1) = ((genes(j) - means2(j, 1))^2) + ((genes(j) - means2(j, 2))^2);
          
        end


        % calculate the sse of those means (sum it in order to look for the
        % minimum value)
        SSE1(i, 1) = ((genes(i) - means(i, 1))^2) + ((genes(i) - means(i, 2))^2);

    end


    % get the minimum sse values
    [minSSE1, index1] = min(SSE1);
    [minSSE2, index2] = min(SSE2);


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
    F1 = MSR1/MSE1;

    % F statistic for two step
    F2 = MSR2 / MSE2;

    % F statistic represents the relative goodness of fit of a one-step 
    % versus a two-step pattern.
    F12 = ((minSSE1 - minSSE2) / (4 - 3)) / (minSSE2 / (n - 4));

    % decide if it is one, two step or other
    if (F1 < 0.05 && (F12 >= 0.05))
        step = "OneStep";

    elseif (F2 < 0.05)
        step = "TwoStep";
    
    else
        step = "Other";
    end

end