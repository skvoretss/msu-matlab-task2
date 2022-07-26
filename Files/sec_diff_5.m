function f = sec_diff_5(x)
f = (200./(x.^3)).*(50.*sin(100./x)./x - cos(100./x));
end