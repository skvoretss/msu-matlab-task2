% опорна€ функци€ равна sup(x*t1 + y*t2) по всем (x, y) из множества
%A, Aeq - матрицы, b, beq - вектора, c, ceq - функции, возвращающие вектор,
%f(x) возвращает скал€р, x0, lb, ub - вектора или матрицы, 
function [func] = supportLebesgue(f, x0, A, b, Aeq, beq, lb, ub)
scalar = @(x, t) x(1)*t(1) + x(2)*t(2);
function [c, ceq] = nlcon(x)
    c = f(x);
    ceq = [];
end
nonlcon = @nlcon;
fmin = @(t) fmincon(@(x) -scalar(x,t), x0, A, b, Aeq, beq, lb, ub, nonlcon);
func = @(t) [fmin(t), scalar(fmin(t),t)]; % тут € помен€ла значени€ теперь сначала идет точка, потом значение функции
end