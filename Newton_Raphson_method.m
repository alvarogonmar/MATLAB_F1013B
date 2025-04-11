function resultado = newtonRaphson_method(f, df, x0, tol)
    max_iter = 1000;   % Límite máximo de iteraciones
    iter = 0;          % Contador de iteraciones
    error = inf;       % Inicializa el error como infinito
    
    while error >= tol && iter < max_iter
        fx = f(x0);
        dfx = df(x0);
        
        if dfx == 0
            disp('La derivada es cero. No se puede continuar.');
            resultado = [NaN, NaN];
            return;
        end
        
        x1 = x0 - fx / dfx;
        error = abs((x1 - x0) / x1);
        
        x0 = x1;
        iter = iter + 1;
    end
    
    resultado = [iter, x0];  % Devuelve número de iteraciones y la raíz
end

% Definir función y derivada
f = @(x) x^3 - x - 2;
df = @(x) 3*x^2 - 1;

% Valores iniciales
x0 = 1.5;
tol = 1e-4;

% Llamada a la función externa
resultado = newtonRaphson_method(f, df, x0, tol);

% Mostrar resultados
fprintf('Número de iteraciones: %d\n', resultado(1));
fprintf('Raíz aproximada: %.6f\n', resultado(2));
