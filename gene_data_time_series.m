% get the path of the csv file
path = "GSE71562_E14R012_raw_counts.xlsx";

% load table
data_matrix = readtable(path);

% number of rows
n = size(data_matrix,1);

% remove first column with names
data_matrix(:,1) = [];

% go throught each row and make binary matrix using BASCA

for i = 1:n
    tic;
    binary_matrix(i,:) = BASC_A(data_matrix{i,:});
    T(i) = toc;

end

fprintf("This is the time it took BASC A: %f\n", sum(T))

% get cluster assignment k_means and k_means_op

for i = 1:n
    tic;
    k_means_binary(i,:) = K_Means(2, data_matrix{i,1:18});
    K(i) = toc;

    %tic;
    %k_means_op_binary(i,:) = K_Means_Op(2, data_matrix{i,1:18});
    %K_Op(i) = toc;
    
end

fprintf("This is the time difference of Kmeans and KmeansOp: %f\n", sum(K))

% compare k_means
%[clusters, count] = groupcounts(result');
%[clusters2, count2] = groupcounts(result2');

%count'
%clusters'

%count2'
%clusters2'

