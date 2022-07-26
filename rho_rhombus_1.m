function [point, val] = rho_rhombus_1(params)
    a = 1; b = 1; x0 = 0; y0 = 0;
    x = [0 0];
    if isfield(params, 'a')
        a = params.a/2; % a - length of horizontal chord
    end
    if isfield(params, 'b')
        b = params.b/2; % b - length of vertical chord
    end
    if isfield(params, 'x0')
       x0 = params.x0;
    end
    if isfield(params, 'y0')
       y0 = params.y0;
    end
    if isfield(params, 't')
       x = params.t;
    end
    val =  max(abs(x(1))*b, abs(x(2))*a) + x(1)*x0 + x(2)* y0;
    point = [x0 + (b*abs(x(1))>=a*abs(x(2)))*b*sign(x(1)), y0 + (b*abs(x(1))<a*abs(x(2)))*a*sign(x(2))];
end