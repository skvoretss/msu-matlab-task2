%��������� �������� � 0 �� [0,pi], �� ��� ������������������ ����������
function f = fPoint(n, x)
f = zeros(size(x));
for i = 1:length(x)
    if (x(i) > 0) && (x(i) < pi/n)
        f(i) = 1;
    end
end
end