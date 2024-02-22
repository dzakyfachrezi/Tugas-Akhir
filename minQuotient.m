function x = minQuotient(m,n)
% Divide m into n numbers with minimum difference

x = zeros(n,1);

for i = 1:n
    x(i) = floor(m/n);
end

for i = 1:mod(m,n)
    x(i) = x(i) + 1;
end

x = sort(x);

end