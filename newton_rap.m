function [r, cont] = newton_rap(f, x0, tol)

syms x;
g = matlabFunction(diff(f, x)); % Derivada de f
counter = 0;

while true
    counter = counter + 1;
    x_new = x0 - f(x0)/g(x0); % Fórmula de Newton-Raphson
    error = abs((x_new - x0) / x_new); % Error relativo
    x0 = x_new;
    if error < tol
        break;
    end
end

r = x0; % Resultado final
cont = counter; % Número de iteraciones

end