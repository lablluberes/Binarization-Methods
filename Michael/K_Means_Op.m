% K-means algorithm optimized?

% This k-means functions has parameters as k (number of clusters)
% and genes (the genes dataset) it returns an array with flags of each
% genes and their assignment in the clusters

% This version calculates the new mean of the centers points 
% everytime a new point is assigned to the cluster. 

% How to run: 
% K_Means_Op(k = num of clusters, genes = dataset)
function result = K_Means_Op(k, genes)

        % Calculates the new mean of the centers points 
        % everytime a new point is assigned to the cluster.

        %tic;

        % get size of the array (columns)
        n = size(genes,2);
    
        % get size of the array (rows)
        m = size(genes,1);
   
    
        % this is the current cluster assignment
        % randomly created with n zeros
        cluster_curr = zeros(1, n);
    
                    %index = randi(n, 1, k);
         index = randperm(n,k);
        % create center_points using rand
         center_points = genes(index);
                     %center_points = genes([1:k], :);
    
                     %center_points = rand(k, m);
    
                     %cluster_curr([1:k]) = [1:k];
        
        center_prev = 0;
    
        % while loop for the algorithm, if the curr and previous center points do
        % not change then stop
        while ~isequal(center_prev, center_points)
    
            % copy the data from the most recent center points to the
            % previous center points
            center_prev = center_points;
            
            % for each data point (gene)
            for i = 1:n
    
                % for each cluster number (k)
                for j = 1:k
    
                    % calculate the distance of cluster k with the gene and
                    % assign it to a variable called distance
                    distance(j) =  abs(center_points(j) - genes(i));
                    %sqrt(sum((center_points(j,:) - genes(i,:)).^2));
    
                end
    
                % find the minimum distance number and index 
                [min_distance, min_index] = min(distance);
    
                % save in the current cluster assigment variable the index of
                % that assigned gene
                cluster_curr(i) = min_index;
    
                % create new center point for another iteration by calculating
                % the mean of the assigned genes if k cluster
                center_points(min_index) = mean(genes(cluster_curr == min_index));
    
            end
    
            %figure;
            %gscatter(genes(:,1), genes(:,2), cluster_curr);
            %hold on
            %scatter(center_points(:,1), center_points(:,2), 100, 'x', 'black', 'linewidth', 3)
    
        end
    
        %figure;
        %gscatter(genes(:,1), genes(:,2), cluster_curr);
        %hold on
        %scatter(center_points(:,1), center_points(:,2), 100, 'x', 'black', 'linewidth', 3)
      
    result = cluster_curr;


end