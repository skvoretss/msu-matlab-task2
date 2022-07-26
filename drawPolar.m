function drawPolar(rho, N, params)
figure;
grid on;
t = params.t;
v = params.c;
t = [t(1)+v(1) t(2)+v(2)];
x = t(1)/sqrt(t(1)^2 + t(2)^2);
y = t(2)/sqrt(t(1)^2 + t(2)^2);

t = [x y];
thi = 2*pi/N;
C = [cos(thi), -sin(thi); sin(thi), cos(thi)];
polar = zeros(N, 2);
params.t = t;
for i = 1:N
    z = rho(params);
    if i == 1
        res01 = z(3);
        v_0x = z(1);
        v_0y = z(2);
        v0x = z(1);
        v0y = z(2);
        if res01 > 0
            polar(i,:) = t./res01;
        elseif res01 < 0
            new_p = params;
            new_p.t = t.*(-1);
            [vec, val] = rho(new_p);
            polar(i,:) = vec.*(-1/val);
        end
    else
        res_1 = z(3);
        v_x = z(1);
        v_y = z(2);
        hold on;
        line([v0x, v_x], [v0y, v_y], 'Color', 'red');
        if res_1 > 0
            polar(i,:) = t./res_1;
        elseif res_1 < 0
            new_p = params;
            new_p.t = t.*(-1);
            z = rho(new_p);
            vec(1) = z(1);
            vec(2) = z(2);
            val = z(3);
            polar(i,:) = vec.*(-1/val);
        end
        v0x = v_x;
        v0y = v_y;
    end
    t = C*t'; %��������b ��� �� ���� thi
    t = t';
    params.t = t;
end
line([v_0x, v_x], [v_0y, v_y], 'Color', 'red'); % ��������� ��������� ����� � ������
k = convhull(polar);
array = polar(k, :);
patch(array(:,1), array(:,2), 'green');
hold off;
end

