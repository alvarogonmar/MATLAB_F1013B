t = linspace(0, 2*pi, 100);
s = linspace(0, 2*pi, 100);

X = [];
Y = [];
Z = [];

for i = 1:length(t)
    for j = 1:length(s)
        X(i, j) = 3*cos(t(i)) + cos(t(i))*cos(s(j));
        Y(i, j) = 3*sin(t(i)) + sin(t(i))*cos(s(j));
        Z(i, j) = sin(s(j));
    end
end

figure;
plot3(X, Y, Z, 'b');
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Curva Parametrizada usando bucles for');
