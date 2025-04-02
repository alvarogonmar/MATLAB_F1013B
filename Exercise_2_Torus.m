% Definición de los parámetros
t = linspace(0, 2*pi, 100); % Valores de t
s = linspace(0, 2*pi, 100); % Valores de s

% Ecuaciones paramétricas con ajuste en Z
X = 3*cos(T) + cos(T).*cos(S);
Y = 3*sin(T) + sin(T).*cos(S);
Z = sin(S); % Escalar para limitar Z entre -0.5 y 0.5

% Graficar la curva parametrizada
figure;
plot3(X(:), Y(:), Z(:), 'b'); % Usamos plot3 para graficar en 3D
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Curva Parametrizada con Z ajustado');


