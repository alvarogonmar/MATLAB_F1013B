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