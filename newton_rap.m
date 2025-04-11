function [r, count] = newton_rap(f, x0, tol)

syms x;
g = matlabFunction(diff(f, x)); % Derivada de f
counter = 0;

while True
    counter = counter + 1;
    % Verificar si f(x0) es igual a 0
    if f(x0) == 0
        r = x0; % Si f(x0) = 0, la raíz es x0
        count = counter;
        return;
    end
    x_new = x0 - (f(x0)/g(x0)); % Fórmula de Newton-Raphson
    error = abs((x_new - x0) / x_new); % Error relativo
    x_new = x0;
    if error < tol
        break;
    end
    x0 = x_new; % Actualizar x0 para la próxima iteración

end

r = x0; % Resultado final
count = counter; % Número de iteraciones