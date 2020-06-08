
function f = func(x)    
    [n,~] = size(x);
    f = 0;
    for i = 1:(n-1)
        for j = (i+1):n
            f = f + 1/norm(x(i,:) - x(j,:));
        end
    end
end