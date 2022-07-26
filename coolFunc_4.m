function f = coolFunc_4(x)
f = (x < 10).*(x.^2) + ((x >= 10) & (x <= 30)).*((x.^2).*(-1))+(x > 30).*(x.^2);
end