%����������������� �������� � 0 �� [0,1], �� ��� ���������� ����������
function f = fRoot(n, x)
f = zeros(size(x));
for i = 1:length(x)
    if (x(i) > 5/n) && (x(i) < 6/n)
        f(i) = 1;
    end
end
end