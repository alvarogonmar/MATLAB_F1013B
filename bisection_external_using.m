clc;
close all;
clear;

f = @(x) sqrt(3*x) - 4;
a = -10;
b = 10;
tol = 0.5;

result = bisection_simple_external(f, a, b, tol);
fprintf('Root found at x = %.4f after %d iterations\n', result(2), result(1));
