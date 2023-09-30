% binarization method based on PPI paper on google drive

function binary = PPI(genes)

    % genes mean
    U = mean(genes);

    % standard deviation
    std_dev = sum((U - genes).^2) / size(genes,2);

    % volatility of gene
    V = (1 / (1 + std_dev)) * 4;

    % threshold 
    G = U + 2 * sqrt(std_dev) * V;

    % binarize based on threshold
    for i = 1:size(genes,2)
        if(genes(i) <= G)
            binary(i) = 0;
        else
            binary(i) = 1;
        end
    end

end