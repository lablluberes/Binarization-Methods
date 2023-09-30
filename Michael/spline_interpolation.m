path = "../HIVIn(Matlab).csv";
path2 = "../LeukB4(Matlab).csv";

data_hiv = readtable(path);
data_leuk = readtable(path2);

data_leuk = rmmissing(data_leuk);

for j = 1:size(data_leuk, 1)

    n_hiv = length(data_leuk{j,:});
    
    x = 1:n_hiv;
    y = data_leuk{j,:};
    
    plot(x, y, 'DisplayName', 'Original', 'LineWidth',1);
    hold on
    title("BASC_A Gene " + j);
    
    for i = [10 20 40 80]%10:10:40
    
        %xx = (n_hiv - 1).*rand(1, i) + 1;
        xx = 1:((n_hiv-1)/(i-1)):n_hiv;
        xx = sort(xx);
    
        yy = spline(x, y, xx);
    
        %thr = mean(K_Means(2, yy));
        %[~, ~, thr] = Stepminer(yy);
        [~, thr] = BASC_A(yy);
    
        yline(thr, '-.', "Thr " + i, 'DisplayName', "Thr " + i);
        plot(xx, yy, '--', 'DisplayName', "Sample " + i );
    end
    
    legend('Location','northwest');
    hold off

    f = gcf;

    %exportgraphics(f,"images/leuk/kmeans/kmeans_gene"+j+".png")
    %exportgraphics(f,"images/leuk/stepminer/stepminer_gene"+j+".png")
    exportgraphics(f,"images/leuk/basc/basc_gene"+j+".png")
end