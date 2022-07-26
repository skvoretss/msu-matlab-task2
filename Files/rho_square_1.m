function  [point, val] = rho_square_1(params)
    r = 1; x0 = 0; y0 = 0; x = [0 0];
    if isfield(params, 'r')
        r = params.r; % r - length of a side
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
    val = r*(abs(x(1)) + abs(x(2)))/2 + x(1)*x0 + x(2)*y0;
    point = [(x(1)>0)*r/2 + (x(1)<=0)*(-r/2)+x0, (x(2)>0)*r/2 + (x(2)<=0)*(-r/2)+y0];
end 