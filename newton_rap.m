function [r, count] = newton_rap(f, x0, tol)

syms x;
g = matlabFunction(diff(f, x));
counter = 0;

while True
    counter = counter + 1;
    x_new = x0 - (f(x0)/g(x0));
    if
        break
    end

end

count = counter