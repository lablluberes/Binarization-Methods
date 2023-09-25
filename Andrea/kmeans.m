%k is centroids
function [cluster, t] = kmeans(x,k)

%length of vector x
n = length(x);
unchanged = 0;


cluster = zeros(1, n);
distances = zeros(1,k);

%index randomization
z = randperm(n,k);

%array of k's
Ks = x(z);



while ~isequal(unchanged, Ks)

    %update value to compare k mean changes
    unchanged = Ks;

    %for each value in x
    for i = 1:n
        %for each center k
        for j = 1:k
            distances(j) = abs(x(i)-Ks(j));
        end

        %assign index in Ks to array of clusters
        %to group them
        [~, index] = min(distances);
        cluster(i) = index;
    end



    %new means
    %for every k find mean of values attributed to that k group
    for i = 1:k
        Ks(i) = mean(x(cluster == i));
    end
    
end

t = mean(Ks)

end