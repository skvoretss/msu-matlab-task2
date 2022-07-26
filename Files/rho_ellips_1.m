function [val, point] = rho_ellips_1(params) %x in R^2
    a = 1; b = 1; x0 = 0; y0 = 0; alpha = 0;
    if isfield(params, 'a')
        a = params.a;
    end
    if isfield(params, 'b')
       b = params.b;
    end
    if isfield(params, 'x0')
       x0 = params.x0;
    end
    if isfield(params, 'y0')
       y0 = params.y0;
    end
    if isfield(params, 'alpha')
       alpha = params.alpha;
    end
    T = [a*a, 0; 0, b*b];
    C = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
    P = C'.*T.*C;
    Pl = P*x';
    f = @(x) sqrt(x(1)*Pl(1) + x(2)*Pl(2)) + x(1)*x0 + x(2)*y0;
    val = f(x);
    point = [x0 y0] + Pl'./sqrt(x(1)*Pl(1) + x(2)*Pl(2));
end