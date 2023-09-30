path = "../HIVIn(Matlab).csv";
path2 = "../LeukB4(Matlab).csv";

data_hiv = readtable(path);
data_leuk = readtable(path2);

data_leuk = rmmissing(data_leuk);

n_hiv = size(data_hiv,1);
%n_leuk = size(data_leuk,1);

%for i = 1:n_hiv
%    hiv_binary(i,:) = BASC_A(data_hiv{i,:});
%end

% problem with dataset there is one gene with nan value
%for i = 1:n_leuk
%    leuk_binary(i,:) = BASC_A(data_leuk{i,:});
%end


%for i = 1:n_hiv
%     hiv_kmeans(i,:) = K_Means(2, data_hiv{i,:});
%end

%for i = 1:n_hiv
%     hiv_kmeans_op(i,:) = K_Means_Op(2, data_hiv{i,:});
%end


%for i = 1:n_hiv
%     [binary, step, ~] = Stepminer(data_hiv{i,:});
%     hiv_step(i,:) = binary;
%end

