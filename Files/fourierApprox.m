function fourierApprox(f, a, b, n, meth)
x = linspace(-1, 1);
x_new = linspace(a, b);
j = 0;

y_f = f(x_new);
%disp(y_f);
fig = figure;

if strcmp(meth, 'Legendre')
    j = 0;
    s = 0.5 * trapz(x, y_f); % i == 0
    for i = 1:n
        f_i = getFuncLgnd(i);
        Y = y_f .* f_i(x);
        c_i = (2*i+1)/2 * trapz(x, Y); 
        s = s + c_i.*f_i(x);
        pause(0.015);
        j = j + 1;
        M(j) = getframe(fig);

        plot(x_new, s, '--', x_new, y_f);
        title(meth);
        xlabel('Abscissa');
        ylabel('Ordinate');
        legend('i-sum', 'f');
    end
elseif strcmp(meth, 'Chebyshev')
    x = linspace(-1+0.0001, 1-0.0001);
    s = 1/pi * trapz(x, y_f./(1 - x.^2).^(0.5)); % i == 0
    j = 0;
    for i = 1:n
        f_i = getFuncCheb(i);
        Y = y_f .* f_i(x) ./ (1 - x.^2).^(0.5);
        c_i = 2/pi * trapz(x, Y);
        s = s + c_i.*f_i(x);
        pause(0.015);
        j = j + 1;
        M(j) = getframe(fig);
        plot(x_new, s, '--', x_new, y_f);
        title(meth);
        xlabel('Abscissa');
        ylabel('Ordinate');
        legend('i-sum', 'f');
    end
elseif strcmp(meth, 'Trigonometry')
    j = 0;
    x = linspace(a, b);
    y_f_n = f(x);
    s = 0.5 * trapz(x, y_f_n) * 1/b; % i == 0
    for i = 1:n
        f_i_cos = @(x) cos(pi*i*x/b);
        f_i_sin = @(x) sin(pi*i*x/b);
        Y_cos = y_f_n .* f_i_cos(x);
        Y_sin = y_f_n .* f_i_sin(x);
        a_i = 1/b * trapz(x, Y_cos);
        b_i = 1/b * trapz(x, Y_sin);
        s = s + a_i .* f_i_cos(x) + b_i .* f_i_sin(x);
        pause(0.015);
        j = j + 1;
        M(j) = getframe(fig);
        plot(x_new, s, '--', x_new, y_f);
        title(meth);
        xlabel('Abscissa');
        ylabel('Ordinate');
        legend('i-sum', 'f');
    end
end
%{
v = VideoWriter(meth);
v.FrameRate = 40;
open(v);
for j = 1:length(M)
    writeVideo(v, M(j));
end
close(v);
%}
disp('NO PROBLEM');
end