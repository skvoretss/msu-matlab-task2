function drawBall(alpha, level, params)
x = linspace(params.ax, params.bx, params.h);
y = linspace(params.ay, params.by, params.h);
z = linspace(params.az, params.bz, params.h);
[X, Y, Z] = meshgrid(x, y, z);
if isinf(alpha)
    V = max(max(abs(X), abs(Y)), abs(Z));
else
    V = abs(X).^(alpha) + abs(Y).^(alpha) + abs(Z).^(alpha);
end
i = isosurface(X,Y,Z,V,level);
if (size(i.vertices) == 0)
    error('SET IS EMPTY');
end

p = patch(i);
xlabel('Abscissa');
ylabel('Ordinate');
zlabel('Applicate');
title('Ball');
grid on;
isonormals(X, Y, Z, V, p);
p.FaceColor = params.faceColour;
p.EdgeColor = params.edgeColour;
daspect([1 1 1])
colormap spring;
view(3); 
axis tight
camlight('left');
%lighting gouraud;


%set(p, 'EdgeColor', 'interp');
disp("NO PROBLEM");
end