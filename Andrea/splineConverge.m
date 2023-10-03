function t = splineConverge(array, n, fun, a)

    pp = spline(1:n, array);
    doesConverge = false;
    counter = 1;
    t1 = -1;
    while ~doesConverge && counter <= 2000

    yy = ppval(pp,1:n/(10*counter):n);

    if fun == "kmeans"
       [~, thresholds(counter)] = kmeans(yy,2);
    elseif fun == "BASC"
       [~, thresholds(counter)] = BASCA(yy);
    elseif fun == "step"
       [~, thresholds(counter)] = stepminer(yy,0.05);
    end
    
    if(length(thresholds) > 1)
        for i = 1:length(thresholds)
            for j = 1:length(thresholds)
                if j ~= i
                    
                    if abs(thresholds(j) - thresholds(i)) <= a
                        doesConverge = true;
                        t = thresholds(j);
                        break
                    end
                end                
            end    
            if doesConverge
                break
            end
        end

    end

    counter = counter + 1;

    end

    counter
end