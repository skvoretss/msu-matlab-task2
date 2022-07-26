function res = getFuncLgnd(n)
if n == 0
    res = @(x) 1;
elseif n == 1
    res = @(x) x;
else
    P1 = getFuncLgnd(n-1);
    P2 = getFuncLgnd(n-2);
    res = @(x) (2*n - 1).*x.*P1(x)./n - (n-1).*P2(x)./n;
end

end