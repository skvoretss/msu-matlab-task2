function drawSetLeb(rho, N, params)
figure;
x_in = zeros(1, N + 1);
y_in = zeros(1, N + 1);
x_out = zeros(1, N + 1);
y_out = zeros(1, N + 1);

t = params.t;
v = params.c;
thi = 2*pi/N;

C = [cos(thi), -sin(thi); sin(thi), cos(thi)];

a1 = 0; b1 = 0; c1 = 0;

a0 = 0; b0 = 0; c0 = 0;

t = [t(1)+v(1) t(2)+v(2)];
x = t(1)/sqrt(t(1)^2 + t(2)^2);
y = t(2)/sqrt(t(1)^2 + t(2)^2);

t = [x y];
params.t = t;

for i = 1:N

    z = rho(params.t);

    x_in(i) = z(1);
    y_in(i) = z(2);
    
    if i == 1 % для пересечения последней с последней касательных
        
        a0 = -t(2);
        b0 = t(1);
        c0 = -t(1)*z(1)-t(2)*z(2);
        
        
        a1 = -t(2);
        b1 = t(1);
        c1 = -t(1)*z(1)-t(2)*z(2);

    else
        a2 = -t(2);
        b2 = t(1);
        c2 = -t(1)*z(1)-t(2)*z(2);
        
        y_out(i-1) = (b1*c2 - b2*c1)/(b2*a1 - a2*b1);
        x_out(i-1) = (-a1*c2 + a2*c1)/(b2*a1 - a2*b1);
        
        a1 = a2;
        b1 = b2;
        c1 = c2;

    end
    t = C*t'; %крутануль еще на угол thi
    t = t';
    params.t = t;
end

%замкнули 
x_in(N+1) = x_in(1);
y_in(N+1) = y_in(1);

% ищем пересечение первой с последней
y_out(N) = (b0*c1 - b1*c0)/(b1*a0 - a1*b0);
x_out(N) = (-a0*c1 + a1*c0)/(b1*a0 - a1*b0);

%замкнули
x_out(N+1) = x_out(1);
y_out(N+1) = y_out(1);

disp('External coors:');
disp(x_out);
disp(y_out);

disp('Internal coors:');
disp(x_in);
disp(y_in);

plot(x_in, y_in, '--', x_out, y_out, '.-');
axis tight;
grid on;
title(func2str(rho));
xlabel('Abscissa');
ylabel('Ordinate');
legend('Internal', 'External');

end