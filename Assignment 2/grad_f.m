
function g_f = grad_f(x)
    [n,d] = size(x);
    n_bar = n*d;
    g_f = zeros(n_bar,1);
    sum_over = 1:n;
    l = 0;
    for i = 1:n
        for j = 1:d
            l = l+1;
            for k = sum_over(sum_over ~= i)  % sum over all except ith
                g_f(l) = g_f(l) -(x(i,j) - x(k,j))/norm(x(i,:) - x(k,:))^3;
            end
        end
    end
end 