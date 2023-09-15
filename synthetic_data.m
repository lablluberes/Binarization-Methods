% number of rows
n = 300;

% number of columns
m = 10;

% data from paper for BASC
genes = [  0.2200    0.2900    0.1000    0.1300    0.8000    0.9000    0.2200    0.8500    0.8100    0.5000];
genes2 = [  0.5000    0.5700    0.1000    0.1300    0.8000    0.9000    0.5000    0.8500    0.8100    0.5800];

% data for cluster
genes3 = [randn(n,1) + 2 randn(n,1) - 2; randn(n, 1) - 2 randn(n,1) + 2];

% data random for basc
genes4 = [ rand(n, m) ];


%%% Print time execution of algorithms BASC A and both versions of K Means

bascA = @() BASC_A(genes);
fprintf("Time of BASC A using genes: %f\n", timeit(bascA))

bascA = @() BASC_A(genes2);
fprintf("Time of Basc A using genes2: %f\n", timeit(bascA))

for i = 1:n
    tic;
    binary_matrix(i,:) = BASC_A(genes4(i,:));
    T(i) = toc;
    
end

fprintf("This is the time it took BASC A on genes4: %f\n", sum(T))

% get cluster assignment k_means and k_means_op

for i = 1:n
    tic;
    k_means_binary(i,:) = K_Means(2, genes4(i,:));
    K(i) = toc;

    tic;
    k_means_op_binary(i,:) = K_Means_Op(2, genes4(i,:));
    K_Op(i) = toc;
    
end

fprintf("This is the time difference of Kmeans and KmeansOp: %f\n", sum(K) - sum(K_Op))




