function convergenceFunc(fn,f,a,b, n, convType)
fig = figure;

if strcmp(convType, 'Pointwise convergence')
    s = '';
end
if strcmp(convType, 'Uniform convergence')
    s = 'sup||fn - f||';
end
if strcmp(convType, 'Root-mean-square convergence')
    s = '||f||_2';
end

m = 100;
x = linspace(a, b, m); 
y = f(x);

j = 0;
y_min = 0;
y_max = 0;
for i = 1:n
    pause(0.015);
    j = j + 1;
    M(j) = getframe(fig);
    y_i = fn(i, x);
    if i == 1
        y_min = min(y);
        y_max = max(y_i);
        if y_min == y_max
            y_max = Inf;
        end
    end
    plot(x, y_i, '--', x, y);
    axis ([a b y_min y_max]);
    if strcmp(convType, 'Uniform convergence')
        title(string(max(abs(y_i - y))));
    elseif strcmp(convType, 'Root-mean-square convergence')
        r = sqrt(trapz(x, abs(y_i - y)));
        title(string(r));
    end
    xlabel('Abscissa');
    ylabel('Ordinate');
    legend(func2str(fn), 'f');
    
end

v = VideoWriter(func2str(fn));
v.FrameRate = 40;
open(v);
for j = 1:length(M)
    writeVideo(v, M(j));
end
close(v);
disp('NO PROBLEM');
