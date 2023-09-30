% This is the code made for the BASC A algorithm

% How to run:
% BASC_A(genes = vector) The vector does not need to be ordered. This is
% done in the code

function [result, thr] = BASC_A(genes)

    originalGenes = genes;
    genes = sort(originalGenes);

    %%% Step 1: Compute a Series of Step Function

    % size of the vector of genes
    n = size(genes,2);

    % initialization C_i_(0) = c_i_N
    % calculate first cost matrix column with no intermidiate break points
    for i = 1:n
        cost_matrix(i, 1) = C_a_b(genes, i, n);
    end


     % Algorithm 1: Calculate optimal step functions
    for j = 1:n-2
        for i = 1:n-j

            % find the min_value and index of the minimum cost of
            % additional intermediate break points
            min_value = Inf;
            min_index = Inf;
            
            for d = i:n-j
    
               curr_value = C_a_b(genes, i, d) + cost_matrix(d+1,j);
          
               if (curr_value < min_value)
                   min_value = curr_value;
                   min_index = d;
               end

            end

            cost_matrix(i, j + 1) = min_value;
            ind_matrix(i, j) = min_index;
            
        end
    end

    cost_matrix;
    ind_matrix;

    %  Algorithm 2: Compute the break points of all optimal step functions

    for j = 1:n-2
        z = j;
        P(1,j) = ind_matrix(1,z);

        if (j > 1)
            z = z - 1;

            for i = 2:j
                P(i,j) = ind_matrix(P(i-1,j)+1, z);
                z = z - 1;
            end
        end
    end

    P;

    %%% Step 2: Find Strongest Discontinuity in Each Step Function

    v = zeros(1, n-2);

    % go through P and find the biggest discontinuity 
    for j = 1:size(P,2)
    
        max_value = -Inf;
        max_index = j;

        % calculate q_score using h:height and e:error
        % the biggest discontinuity then is saved in v; 
        for i = 1:j

            h = determine_h(P, i, j, genes);
  
            z = ((genes(P(i, j))) + genes((P(i, j)+1))) / 2;
        
            e = sum((genes(1:n) - z).^2);

            q_score = h / e;

            if(q_score > max_value)

                max_value = q_score;
                max_index = i;

            end

        end

        v(j) = P(max_index, j);

    end

    v;
    
    %%% Step 3: Estimate Location and Variation of the Strongest Discontinuities

    % calculate the threshold using median of v
    thr = (genes(round(median(v))+1) + genes(round(median(v)))) / 2;

    % for each gene assign 0 or 1 based on the threshold
    for i = 1:n
        if(originalGenes(i) <= thr)
            u_binary(i) = 0;
        else
            u_binary(i) = 1;
        end
    end

    % length of v
    n = length(v);

    % normalized vector v
    v0 = v(1:n-2) / (n - 1);

    % average deviation
    AD_v = (1/n) * abs(sum(v(1:n-2) - median(v)));

    t = 0.15;

    % test statistics
    t_0 = t - AD_v;
    
    % return the binary vector
    result = u_binary;
    

end

% function to calculate cost C_a_b
function quadratic_distance = C_a_b(genes, a, b)

    mean = Y_a_b(genes, a, b);
   
    quadratic_distance = sum( (genes(a:b) - mean ).^2 );

end

% function to calculate mean Y_a_b
function mean = Y_a_b(genes, a, b)
        
    mean = sum(genes(a:b) / (b - a + 1) );
    
end

% function to determine the height
function h = determine_h(P, i, j, genes)

    N = size(genes,2);

    if (i == 1 && j > 1)
        h = Y_a_b(genes, P(i,j) + 1, P(i+1,j)) - Y_a_b(genes, 1, P(i,j));

    elseif (i == j && j > 1)
        h = Y_a_b(genes, P(i,j) + 1, N) - Y_a_b(genes, P(i-1,j) + 1, P(i,j));

    elseif (i == 1 && j == 1)
        h = Y_a_b(genes, P(i,j) + 1, N) - Y_a_b(genes, 1, P(i,j));

    else
        h = Y_a_b(genes, P(i,j) + 1, P(i+1,j)) - Y_a_b(genes, P(i-1,j) + 1, P(i,j));
    end
end











