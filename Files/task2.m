%#ok<*NOPTS>
%% Paragraphs 1-2
clc;
a = 1;
h = 0.005;
b = 50;
d_h = 100;
xx = a:h:b;
x = xx(1:d_h:end);
%f = @coolFunc_1;
compareInterp(x, xx, @coolFunc_1);
%f = @coolFunc_2;
compareInterp(x, xx, @coolFunc_2);
%f = @coolFunc_3;
compareInterp(x, xx, @coolFunc_3);
%f = @coolFunc_4;
compareInterp(x, xx, @coolFunc_4);
%f = @coolFunc_5;
compareInterp(x, xx, @coolFunc_5);
%f = @coolFunc_6;
compareInterp(x, xx, @coolFunc_6);
a = 1;
h = 0.005;
b = 5;
d_h = 3;
xx = a+h:h:b;
x = xx(1:d_h:end);
compareInterp(x, xx, @coolFunc_2);
clear;


%% Paragraph 3
clc;
a = 1;
h = 0.05;
b = 5;
d_h = 2;
xx = a+h:h:b;
figure;
Y = zeros(1, length(xx) - 1);
XX = zeros(1, length(xx) - 1);
j = 1;
for i = a+h:h:b
    t = i-h:10:i;
    y_f = abs(sec_diff_2(t));
    y = max(y_f);
    x = find(y_f >= y, 1);
    Y(j) = y*h*h/2;
    XX(j) = t(x);
    j = j + 1;
end
f_g = @coolFunc_2;

X = XX(1:d_h:end);
RES = real_pog(X, XX, f_g);
disp(Y);
disp(RES);
plot(XX, Y, '--', XX, RES, '-.');
title('Huge Error Rate');
xlabel('Abscissa');
ylabel('Ordinate' );
legend('Predicted', 'Real');
%%{
figure;
j = 1;
for i = a+h:h:b
    t = i-h:10:i;
    y_f = abs(sec_diff_6(t));
    y = max(y_f);
    x = find(y_f >= y, 1);
    Y(j) = y*h*h/2;
    XX(j) = t(x);
    j = j + 1;
end
g = @coolFunc_6;
X = XX(a:d_h:end);
RES = real_pog(X, XX, g);
plot(XX, Y, '--', XX, RES, '-.');
title('Small Error Rate');
xlabel('Abscissa');
ylabel('Ordinate' );
legend('Predicted', 'Real');
clear;
%}


%% Paragraph 4
clc;
convType1 = 'Pointwise convergence';
convType2 = 'Uniform convergence';
convType3 = 'Root-mean-square convergence';
a = 1;
b = 10;
n = 100;

%Pointwise convergence
convergenceFunc(@f_n, @f, a, b, n, convType1); 

%Uniform convergence
convergenceFunc(@fUni, @f, a, b, n, convType2); 

a = 0;
b = 0.99;
%Pointwise convergence and no Uniform convergence
convergenceFunc(@fPnu, @f, a, b, n, convType1);  

%Root-mean-square convergence and no Pointwise convergence
convergenceFunc(@fRoot, @f, a, b, n, convType3); 

%% Paragraph 5 - Legendre's system has been reliased 
clc;
fourierApprox(@(x) exp(x), -10, 3, 10, 'Legendre')
fourierApprox(@(x) x, -5, 5, 40, 'Chebyshev')
fourierApprox(@(x) sin(x), -4, 6, 20, 'Trigonometry') % смотри теорему Дирихле 
clear;

%% Paragraph 6
clc;
% Part 1

g = @(x) sin(x);
a = -5;
b = 5;
h = 1000;

%{
g = @(x) (1 - x.^2) .* exp((-1/2)*x.^2);
a = -5;
b = 1;
h = 600;
%}
x = linspace(a, b, h);
% Part 2
y = g(x);
t = 0;
for i = 2:length(y)-1
    if (y(i) < y(i+1)) && (y(i) < y(i-1))
        t = t + 1;
    end
end
j = 1;
y_local_min = zeros(1, t);
x_local_min = zeros(1, t);
index_local_min = zeros(1, t);
for i = 2:length(y)-1
    if (y(i) < y(i+1)) && (y(i) < y(i-1))
        y_local_min(j) = y(i);
        x_local_min(j) = x(i);
        index_local_min(j) = i;
        j = j + 1;
    end
end
y_max = max(y);
x_max = find(y == y_max); % индексы всех максимумов
indx_max = x_max(1); % индекс 1-го максимума
plot(x, y, x_local_min, y_local_min, 'r*', x(indx_max), y_max, 'g*'); % отметили минимумы и максимум и нарисовали график

hold on;
x_indx_min_1 = find(x_local_min > x(indx_max)); % смотрим минимумы справа от максимума
x_indx_min_2 = find(x_local_min < x(indx_max)); % смотрим минимумы слева от максимума

delta_1 = -1;
if ~isempty(x_indx_min_1)
    delta_1 = abs(x_local_min(x_indx_min_1(1)) - x(indx_max));
    disp("x_indx_min_1 is not empty");
end

delta_2 = -1;
if ~isempty(x_indx_min_2)
    delta_2 = abs(x_local_min(x_indx_min_2(1)) - x(indx_max));
    disp("x_indx_min_2 is not empty");
end

if delta_1 ~= -1
    if delta_2 ~= -1
        if delta_1 <= delta_2
            x_fly = x(indx_max:index_local_min(x_indx_min_1(1)));
            y_fly = y(indx_max:index_local_min(x_indx_min_1(1)));
            comet(x_fly, y_fly);
        else
            x_fly = x(index_local_min(x_indx_min_2(1)):indx_max);
            y_fly = y(index_local_min(x_indx_min_2(1)):indx_max);
            comet(flip(x_fly), flip(y_fly));
        end
    else
        x_fly = x(indx_max:index_local_min(x_indx_min_1(1)));
        y_fly = y(indx_max:index_local_min(x_indx_min_1(1)));
        comet(x_fly, y_fly);
    end
else
    if delta_2 == -1
        disp("NO MINIMUM");
    else
        x_fly = x(index_local_min(x_indx_min_2(1)):indx_max);
        y_fly = y(index_local_min(x_indx_min_2(1)):indx_max);
        comet(flip(x_fly), flip(y_fly));
    end
end
hold off;
clear;

%% Paragraph 11 - 12
clc;

a = -10;
h = 1000;
b = 10;
n = 100;
x = linspace(a, b, h);
y = linspace(a, b, h);
[X, Y] = meshgrid(x, y);
Z = sin(X) + cos(Y);

contour(X, Y, Z, 10);

fig = figure;

z_max = Z(2:h-1, 2:h-1) > Z(2:h-1, 1:h-2); % красный квадрат | <- 
z_max = z_max & (Z(2:h-1, 2:h-1) > Z(1:h-2, 2:h-1)); % зеленый квадрат | ? 
z_max = z_max & (Z(2:h-1, 2:h-1) > Z(2:h-1, 3:h)); % синий квадрат | ->
z_max = z_max & (Z(2:h-1, 2:h-1) > Z(3:h, 2:h-1)); % серый квадрат | ?
[i_max, j_max] = find(z_max);

z_min = Z(2:h-1, 2:h-1) < Z(2:h-1, 1:h-2); % красный квадрат
z_min = z_min & (Z(2:h-1, 2:h-1) < Z(1:h-2, 2:h-1)); % зеленый квадрат
z_min = z_min & (Z(2:h-1, 2:h-1) < Z(2:h-1, 3:h)); % синий квадрат
z_min = z_min & (Z(2:h-1, 2:h-1) < Z(3:h, 2:h-1)); % серый квадрат
[i_min, j_min] = find(z_min);

% +1, т.к. мы изначально считали без краёв и нужно вернуть порядок
i_max = i_max + 1; 
j_max = j_max + 1;
i_min = i_min + 1;
j_min = j_min + 1;

% i_max и тд - i-е индексы максимума
% X(i_max, j_max) - матрица, где на диагонали стоят нужные нам значения
% diag(X(i_max, j_max)) - искомый вектор

s = surf(X, Y, Z, 'EdgeColor', 'none');
colormap summer;
hold on;
p_max = plot3(diag(X(i_max, j_max)), diag(Y(i_max, j_max)), diag(Z(i_max, j_max)), 'black *');
p_min = plot3(diag(X(i_min, j_min)), diag(Y(i_min, j_min)), diag(Z(i_min, j_min)), 'r*');
hold off;

j = 0;
delta = 0.04;
for i = 1:n
    %pause(0.01);
    j = j + 1;
    M(j) = getframe(fig);
    %M(j) = getframe();
    s.ZData = sin(X + delta) + cos(Y - delta);
    z_max = s.ZData(2:h-1, 2:h-1) > s.ZData(2:h-1, 1:h-2); % красный квадрат | <- 
    z_max = z_max & (s.ZData(2:h-1, 2:h-1) > s.ZData(1:h-2, 2:h-1)); % зеленый квадрат | ? 
    z_max = z_max & (s.ZData(2:h-1, 2:h-1) > s.ZData(2:h-1, 3:h)); % синий квадрат | ->
    z_max = z_max & (s.ZData(2:h-1, 2:h-1) > s.ZData(3:h, 2:h-1)); % серый квадрат | ?
    [i_max, j_max] = find(z_max);

    z_min = s.ZData(2:h-1, 2:h-1) < s.ZData(2:h-1, 1:h-2); % красный квадрат
    z_min = z_min & (s.ZData(2:h-1, 2:h-1) < s.ZData(1:h-2, 2:h-1)); % зеленый квадрат
    z_min = z_min & (s.ZData(2:h-1, 2:h-1) < s.ZData(2:h-1, 3:h)); % синий квадрат
    z_min = z_min & (s.ZData(2:h-1, 2:h-1) < s.ZData(3:h, 2:h-1)); % серый квадрат
    [i_min, j_min] = find(z_min);

    % +1, т.к. мы изначально считали без краёв и нужно вернуть порядок
    i_max = i_max + 1; 
    j_max = j_max + 1;
    i_min = i_min + 1;
    j_min = j_min + 1;
    
    p_max.XData = diag(s.XData(i_max, j_max));
    p_max.YData = diag(s.YData(i_max, j_max));
    p_max.ZData = diag(s.ZData(i_max, j_max));
    
    p_min.XData = diag(s.XData(i_min, j_min));
    p_min.YData = diag(s.YData(i_min, j_min));
    p_min.ZData = diag(s.ZData(i_min, j_min));
    drawnow;
    delta = delta + 0.04;
    title('Animation');
end

movie(M);

v = VideoWriter('Animation_avi');
v.FrameRate = 40;
open(v);
for j = 1:length(M)
    writeVideo(v, M(j));
end
close(v);


t = matfile('Animation_mat.mat','Writable',true);
t.M = M;

disp('NO PROBLEM');

clear;


%% Paragraph 14
clc;

params.ax = -5;
params.bx = 5; 
params.ay = -5;
params.by = 5;
params.az = -5;
params.bz = 5;
params.h = 100;
params.faceColour = 'yellow';
params.edgeColour = 'none';
alpha = 2;
level = 4;

drawBall(alpha, level, params)
clear;

%% Paragraph 15
clc;
alphas = [4, 2, 1, 0.5, Inf];
colors = ['r', 'g', 'y', 'm', 'b'];
edges = [1, 0, 0, 1, 0];

drawManyBalls(alphas, colors, edges)
clear;


%% Paragraph 7
clc;

getEqual(@(t) sin(t), @(t) cos(t), 0, 10, 10)
getEqual(@(t) sin(3*t), @(t) cos(2*t), 0, 10, 10)
getEqual(@(t) sin(t), @(t) cos(2*t), 0, 10, 10)

%% Paragraph 13
clc;


points = [-10 5; 20 -2; 5 -10; -5 7];
viewPossible(points, 1, 0.2);

points = [-10 5; 20 -2; 5 -10];
viewPossible(points, 0.7, 0.2);
clear;


%% Paragraph 9
clc;

g = @(x) x(1)^2 + x(2)^2 - 1;
x0 = [0, 0]; %начальная точка в середине области
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-3, -3];
ub = [3, 3];
rho = supportLebesgue(g, x0, A, b, Aeq, beq, lb, ub);

t = [13 1];
disp(rho(t));

clear;
%% Paragraph 8
clc;

params.a = 5;
params.b = 4;
params.t = [1 0];
params.c = [0 0];

N = 16;
drawSet(@rho_ellipse, N, params);

drawSet(@rho_rhomb, N, params);

params.a = 2;
params.b = 2;

drawSet(@rho_square, N, params);

g = @(x) x(1)^2 + x(2)^2 - 1;

x0 = [0, 0]; %начальная точка в середине области
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-3, -3];
ub = [3, 3];
rho = supportLebesgue(g, x0, A, b, Aeq, beq, lb, ub);
drawSetLeb(rho, N, params);

N = 100;
params.a = 5;
params.b = 4;
drawSet(@rho_ellipse, N, params);

drawSet(@rho_rhomb, N, params);

params.a = 2;
params.b = 2;

drawSet(@rho_square, N, params);
drawSetLeb(rho, N, params);
%% Paragraph 10
clc;
params.a = 5;
params.b = 4;
params.t = [3 5];
params.c = [0 0];

N = 16;

drawPolar(@rho_ellipse, N, params);
drawPolar(@rho_rhomb, N, params);
params.c = [3 3];

drawPolar(@rho_ellipse, N, params);
drawPolar(@rho_rhomb, N, params);

params.a = 2;
params.b = 2;
params.c = [0 0];
drawPolar(@rho_square, N, params);

params.c = [-3 4];
drawPolar(@rho_square, N, params);