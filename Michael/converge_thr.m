
data = readtable("HIVIn(Matlab).csv");
data2 = readtable("LeukB4(Matlab).csv");
data2 = rmmissing(data2);

% hiv genes
up = data{36,:};
down = data{30,:};
up_down = data{28,:};
down_up = data{5,:};

% leuk genes
leuk_up = data2{10,:};
leuk_down = data2{4,:};
leuk_up_down = data2{380,:};
leuk_down_up = data2{205,:};

n_hiv = length(leuk_up);

converge = false;

samples = 10;
index = 1;

tolerance = 0.00001;


while(~converge)

    x = 1:n_hiv;

    xx = 1:((n_hiv-1)/(samples-1)):n_hiv;
    xx = sort(xx);
    
    yy = spline(x, leuk_up, xx);

    %thr = mean(K_Means(2, yy));
    %[~, ~, thr] = Stepminer(yy);
    [~, thr] = BASC_A(yy);

    threshold(index) = thr;
    index = index + 1;

    n_thr = length(threshold);

    samples = samples + 10;
    %distance = zeros(threshold, threshold);

    for i = 1:n_thr
        for j = 1:n_thr
            if(i ~= j)
                difference = abs(threshold(i) - threshold(j));

                if(difference <= tolerance)

                    t1 = threshold(i);
                    t2 = threshold(j);
                    converge = true;

                    break;
                end
            end
        end
        if(converge)
            break;
        end
    end
    
end

t1
t2
samples

