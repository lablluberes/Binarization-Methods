function [u, t] = BASCA(x)


n = length(x);

sorted = sort(x);

%euclidean distance


%break point indexes

%this is a matrix nxn-1

C = zeros(n,n-1);
ind = zeros(n-1,n-2);
P = zeros(n-2,n-2);

%initialize cost matrix first column

for i = 1:n
    C(i,1) = C_ab(sorted,i,n);
end

%iteration

for j = 1:n-2

    for i = 1:n-j

        cost = Inf;
        %d thru N-j to find mins
        for d = i:n-j
            value = (C_ab(sorted,i,d) + C(d+1,j));

            if value < cost
                cost = value;
                index = d;
            end
        end
        C(i,j+1) = cost;
        %j instead of j+1 here because later ind(i,1) is accessed
        ind(i,j) = index;
    end

end




%breakpoints of stepfunctions


for j = 1:n-2
    z = j
    P(1,j) = ind(1,z);

    if j > 1
        z = z - 1;
        for i = 2:j
            P(i,j) = ind(P(i-1,j)+1,z);
            z = z-1;
        end
    end
end




%discontinuity search

v = zeros(1,n-2)

for j = 1:n-2         

            %discontinuity h

            maxQ = -Inf;
            index = j;

            for i = 1:j

                if i == 1 && j > 1
                    h = Y_ab(sorted,P(i,j)+1,P(i+1,j)) - Y_ab(sorted,1,P(i,j));
                elseif i == j && j>1
                    h = Y_ab(sorted,P(i,j)+1,n) - Y_ab(sorted,P(i-1,j)+1,P(i,j));
                elseif i == j && j == 1
                    h = Y_ab(sorted,P(i,j)+1,n) - Y_ab(sorted,1,P(i,j));
                else
                    h = Y_ab(sorted,P(i,j)+1,P(i+1,j)) - Y_ab(sorted,P(i-1,j)+1,P(i,j));
                end


                z = (sorted(P(i,j)) + sorted(P(i,j)+1))/2;

                e = sum((sorted(1:n)-z).^2);

                q = h/e;

                if q > maxQ

                   maxQ = q;
                   index = i;
                end
            end

            v(j) = P(index,j);
end




%location of strongest discontinuities, threshold

%in case result is a decimal
vmed = median(v);
vmed = round(vmed);
u = zeros(1,n);

t = (sorted(vmed+1) + sorted(vmed))/2;


u = x > t;


end


function val = C_ab(x,a,b)
    

    Y = Y_ab(x,a,b);

    val = sum((x(a:b)-Y).^2);

end

function val = Y_ab(x,a,b)

    val = sum((x(a:b))/(b-a+1));

end