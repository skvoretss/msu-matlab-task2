function res = real_pog_2(x, xx, f)
y_1 = f(xx);
y_2 = f(x);
r = interp1(x, y_2, xx, 'linear');
end
