function res = real_pog(x, xx, f)
y = f(xx);
z = f(x);
r1 = interp1(x, z, xx);
disp(y);
disp(r1);
disp(xx);
res = abs(y - r1);
end