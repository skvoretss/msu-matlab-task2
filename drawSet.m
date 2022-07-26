% ошибки:
%1. плохое хранение предыдущего вектора через a0, b0, c0
%2. ошибки в уравнениях

% внесенные изменения: 
% 1) отправка не случайного вектора в rho а движением на угол thi через
% linspace и [cos(thi) sin(thi)]
% 2) хранение пред.значений через вектора t1 и t2, а не через a0, b0 и т.п
% 3) исправлены коэффициенты в уравнениях

function  drawSet(rho, N, params)
    figure;

    thi = linspace(0, 2*pi, N);
    
    x_in = zeros(1, N);            
    y_in = zeros(1, N);  
    
    x_out = zeros(1, N + 1);
    y_out = zeros(1, N + 1);
    
    t2 = [0 0]; 
   
    for i = 1:N
        t1 = t2;              
        t2 = [cos(thi(i)) sin(thi(i))]; 
        params.t = t2;
        z = rho(params);      
        
        x_in(i) = z(1);
        y_in(i) = z(2);

        if (i == 1)
            c2 = -t2(1)*x_in(i) -t2(2)*y_in(i);
            
        else
            c1 = c2;                
            c2 = -t2(1)*x_in(i) - t2(2)*y_in(i);
            
            y_out(i-1) =(t1(1)*c2 - t2(1)*c1)/(-t1(1)*t2(2) + t1(2)*t2(1));
            x_out(i-1) =(t1(2)*c2 - t2(2)*c1)/(-t1(1)*t2(2) + t1(2)*t2(1));
            
        end
    end

    t1 = t2;
    t2 = [cos(thi(1)) sin(thi(1))];
    c1 = c2;
    c2 = -t2(1)*x_in(1) - t2(2)*y_in(1);

    y_out(N) = y_out(1);
    x_out(N) = x_out(1);
    
    %замкнули
    x_out(N+1) = x_out(1);
    y_out(N+1) = y_out(1);

    plot(x_in, y_in, '--', x_out, y_out, '.-');
    axis tight;
    grid on;
    title(func2str(rho));
    xlabel('Abscissa');
    ylabel('Ordinate');
    legend('Internal', 'External');
end