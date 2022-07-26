%%{
function [res] = getEqual(f, g, t0, t1, N)
figure;
if (t1 <= t0)
    error('T1 <= T0');
end
h = (t1-t0)/(N-1);
t = t0:h:t1; % среднее разбиение
%fgt - матрица точек p
fgt(1:N, 1) = f(t);
fgt(1:N, 2) = g(t);
%all_dist - все возможные расстояния
all_dist = pdist2(fgt, fgt); 
%dist - расстояния от 1-го до 2-го, от 2-го до 3-го и т.д. - взяли верхнюю
%диагональ
dist = diag(all_dist,1); 

s = sum(dist)/(N-1);
new_coor = t;
i = 1;
j = 4*N;
disp("Old dist is:");
disp(dist);
while (max(dist) - min(dist)) > 0.01
    w_j = h*N/j;
    if dist(i) < s
        if i == N-1
            if ((new_coor(i) - w_j) < new_coor(i+1)) && ((new_coor(i) - w_j) > new_coor(i-1))
                new_coor(i) = new_coor(i) - w_j;
                fgt(1:N, 1) = f(new_coor);
                fgt(1:N, 2) = g(new_coor);
                all_dist = pdist2(fgt, fgt); 
                dist = diag(all_dist,1);
                s = sum(dist)/(N-1);
            end
        else
            if ((new_coor(i+1) + w_j) < new_coor(i+2)) && ((new_coor(i+1) + w_j) > new_coor(i))
                new_coor(i+1) = new_coor(i+1) + w_j;
                fgt(1:N, 1) = f(new_coor);
                fgt(1:N, 2) = g(new_coor);
                all_dist = pdist2(fgt, fgt); 
                dist = diag(all_dist,1);
                s = sum(dist)/(N-1);
            end
        end
    elseif dist(i) > s
        if i == N-1
            if ((new_coor(i) + w_j) < new_coor(i+1)) && ((new_coor(i) + w_j) > new_coor(i-1))
                new_coor(i) = new_coor(i) + w_j;
                fgt(1:N, 1) = f(new_coor);
                fgt(1:N, 2) = g(new_coor);
                all_dist = pdist2(fgt, fgt); 
                dist = diag(all_dist,1);
                s = sum(dist)/(N-1);
            end
        else
            if ((new_coor(i+1) - w_j) < new_coor(i+2)) && ((new_coor(i+1) - w_j) > new_coor(i))
                new_coor(i+1) = new_coor(i+1) - w_j;
                fgt(1:N, 1) = f(new_coor);
                fgt(1:N, 2) = g(new_coor);
                all_dist = pdist2(fgt, fgt); 
                dist = diag(all_dist,1);
                s = sum(dist)/(N-1);
            end
        end
    end
    if i == N-1
        i = 0;
    end

    i = i + 1;
    j = j + 1;
end
res(1, 1:N) = f(new_coor);
res(2, 1:N) = g(new_coor);
disp("Old coors are:");
disp(t);
t = linspace(t0, t1, 1000);
disp("New dist is:");
disp(dist);
disp("New coors are:");
disp(new_coor);
plot(f(new_coor), g(new_coor), '-o', f(t), g(t)); 
xlabel('f(t) function');
ylabel('g(t) function');

axis equal;
end
