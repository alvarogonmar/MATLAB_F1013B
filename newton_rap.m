function [r, count] = newton_rap(f, x0, tol)

syms x;
g = matlabFunction(diff(f, x)); % Derivada de f
counter = 0;

while True
    counter = counter + 1;
    x_new = x0 - (f(x0)/g(x0)); % Fórmula de Newton-Raphson
    error = abs((x_new - x0) / x_new); % Error relativo
    if error < tol
        break;
    end
    x0 = x_new; % Actualizar x0 para la próxima iteración

end

r = x_new; % Resultado final
count = counter; % Número de iteraciones