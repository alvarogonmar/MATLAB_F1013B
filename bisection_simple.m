clc;
close all;
clear;

f = @(x) sqrt(3*x) - 4;

a = -10;
b = 10;

tol = 0.5;

if f(a) * f(b) >= 0
    error('Function must have opposite signs at a and b');
end

counter = 0;
while true
    m = (a + b)/2; %midpoint
    if f(m) == 0
        break;
    elseif f(a) * f(m) > 0
        a = m;
    else
        b = m;
    end
    m = (a+b)/2;
    root = m;
    error = abs(b-a)/2;
    counter = counter+1;
    if error < tol
        break;
    end
end

% root = (a+b)/2
% fprintf('Root found at x = %.10f after %d iterations\n"
disp(root)