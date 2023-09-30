function binary = PPI(genes)

    U = mean(genes);

    std_dev = sum((U - genes).^2) / size(genes,2);

    V = (1 / (1 + std_dev)) * 4;

    G = U + 2 * sqrt(std_dev) * V;

    for i = 1:size(genes,2)
        if(genes(i) <= G)
            binary(i) = 0;
        else
            binary(i) = 1;
        end
    end

end