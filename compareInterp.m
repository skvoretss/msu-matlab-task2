function compareInterp(x, xx, f)
y = f(xx);
z = f(x);
figure;
r1 = interp1(x, z, xx, 'linear');
subplot(2,2,1);
plot(xx, y, '--', xx, r1, '-.');
title('Interpolation: linear');
xlabel('Abscissa');
ylabel(func2str(f));
legend('Original', 'Interpolation');

r2 = interp1(x, z, xx, 'nearest');
subplot(2,2,2);
plot(xx, y, '--', xx, r2, '-.');
title('Interpolation: nearest');
xlabel('Abscissa');
ylabel(func2str(f));
legend('Original', 'Interpolation');

r3 = interp1(x, z, xx, 'pchip');
subplot(2,2,3);
plot(xx, y, '--', xx, r3, '-.');
title('Interpolation: pchip');
xlabel('Abscissa');
ylabel(func2str(f));
legend('Original', 'Interpolation');

r4 = interp1(x, z, xx, 'spline');
subplot(2,2,4);
plot(xx, y, '--', xx, r4, '-.');
title('Interpolation: spline');
xlabel('Abscissa');
ylabel(func2str(f));
legend('Original', 'Interpolation');

end