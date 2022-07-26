function viewPossible(points, V, L)
figure;
x = -30:1:30;
y = -30:1:30;
N = length(x);
[X, Y] = meshgrid(x, y);
Z = zeros(N, N);
t = ones(N, 2);
t(:, 2) = y;
for i = 1:N
    t(:, 1) = ones(N, 1);
    t(:, 1) = t(:, 1) .* x(i);
    % t := x_i со всеми y
    s = pdist2(t, points);
    s = s + 1;
    s = V./s;
    d = sum(s, 2);
    Z(:, i) = d;
end
m = contourf(X, Y, Z, [L L], '-.black');
hold on;
plot(points(:, 1), points(:,2), '*r');
xlabel('Abscissa');
ylabel('Ordinate');
hold off;
if m(2) >= (size(m,2)-1)
    disp("Simply connected");
else
    disp("Not simply connected");
end