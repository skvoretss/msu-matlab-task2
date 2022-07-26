function drawManyBalls(alphas, colors, edges)
x = linspace(-3, 3, 15);
y = linspace(-3, 3, 15);
z = linspace(-3, 3, 15);
[X, Y, Z] = meshgrid(x, y, z);
if ((length(alphas) ~= length(colors)) || (length(alphas) ~= length(edges)))
    error('Size is not equal');
end

for i = 1:length(alphas)
    figure;
    if isinf(alphas(i))
        V = max(max(abs(X), abs(Y)), abs(Z));
    else
        V = abs(X).^(alphas(i)) + abs(Y).^(alphas(i)) + abs(Z).^(alphas(i));
    end
    
    j = isosurface(X,Y,Z,V,2);
    
    if (size(j.vertices) == 0)
        error('SET IS EMPTY');
    end
    p = patch(j);
    title('Metric ' + string(alphas(i)));
    xlabel('Abscissa');
    ylabel('Ordinate');
	zlabel('Applicate');
    title('Ball');
    grid on;
    isonormals(X, Y, Z, V, p);
    p.FaceColor = colors(i);
    if (edges(i) == 0)
        p.EdgeColor = 'none';
    end
    daspect([1 1 1])
    view(3); 
    axis tight
    camlight('left');
    %lighting gouraud;
end
end