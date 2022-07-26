% ������� ������� ����� sup(x*t1 + y*t2) �� ���� (x, y) �� ���������
%A, Aeq - �������, b, beq - �������, c, ceq - �������, ������������ ������,
%f(x) ���������� ������, x0, lb, ub - ������� ��� �������, 
function [func] = supportLebesgue(f, x0, A, b, Aeq, beq, lb, ub)
scalar = @(x, t) x(1)*t(1) + x(2)*t(2);
function [c, ceq] = nlcon(x)
    c = f(x);
    ceq = [];
end
nonlcon = @nlcon;
fmin = @(t) fmincon(@(x) -scalar(x,t), x0, A, b, Aeq, beq, lb, ub, nonlcon);
func = @(t) [fmin(t), scalar(fmin(t),t)]; % ��� � �������� �������� ������ ������� ���� �����, ����� �������� �������
end