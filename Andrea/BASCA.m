function u = BASCA(x)


n = length(x)



%euclidean distance


%break point indexes
%to do
breakpointIndex = []

%this is a matrix nxn-1

C = zeros(n,n-1)
ind = zeros(n,n-2)
P = zeros(n,n-2)

%initialize cost matrix first column

for i = 1:n
    C(i,1) = C_ab(x,i,n)
end

%iteration

for j = 1:n-2

    for i = 1:n-j

        %d thru N-j to find mins
        costsearch = zeros(1,n-j-i)
        for d = i:n-j
            costsearch(d-i+1) = (C_ab(x,i,d) + C(d+1,j))
        end
        [cost, index] =  min(costsearch)
        C(i,j+1) = cost
        ind(i,j) = index
    end

end




%breakpoints of stepfunctions


for j = 1:n-2
    z = j
    P(1,j) = ind(1,z);

    if j > 1
        z = z - 1
        for i = 2:j
            P(i,j) = ind(P(i-1,j)+1,z)
            z = z-1
        end
    end
end




%discontinuity search

v = zeros(1,n-2)
qVec = zeros(1,n-2)

for j = 1:n-2         

            %discontinuity h

            for i = 1:j

                if i == 1 & j > 1
                    h = Y_ab(x,P(i,j)+1,P(i+1,j)) - Y_ab(x,1,P(i,j))
                elseif i == j & j>1
                    h = Y_ab(x,P(i,j)+1,n) - Y_ab(x,P(i-1,j)+1,P(i,j))
                elseif i == j & j == 1
                    h = Y_ab(x,P(i,j)+1,n) - Y_ab(x,1,P(i,j))
                else
                    h = Y_ab(x,P(i,j)+1,P(i+1,j)) - Y_ab(x,P(i-1,j)+1,P(i,j))
                end


                z = (x(P(i,j)) + x(P(i,j)+1))/2

                e = sum((x(i:n)-z).^2)

                q = h/e

                qVec(i) = q
            end

            [val, index] = max(qVec)
            v(j) = index
end




%location of strongest discontinuities, threshold

%in case result is a decimal
vmed = median(v)
vmed = round(vmed)
u = zeros(1,n)

t = (x(vmed+1) + x(vmed))/2


for i = 1:n

    if x(i) <= t
       u(i) = 0
    else
        u(i) = 1
    end
    
end


end


function sum = C_ab(x,a,b)
    
    sum = 0

    Y = Y_ab(x,a,b)

    for i = a:b   
        sum = sum + (x(i) - Y)^2
    end

    sum

end

function val = Y_ab(x,a,b)

    val = sum(x(a:b))/(b-a+1)

end