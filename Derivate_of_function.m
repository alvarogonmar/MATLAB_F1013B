%% How to calculate a derivative in Matlab

clear;
clc;
close all;

% define an anonymous function
f = @(x) x.^2 - 3*x - 4;

% First create a symbolic variable (since we are deriving with respect x) 
syms x;

% with the internal MATLAB function we calculate the derivative with respect x
derivative = diff(f, x);

% we convert the previous result to an anonymous function so we can evaluate it
g = matlabFunction (derivative);

% print some results
disp(['f(5): ' num2str(f(5)) ', g(5): ' num2str(g(5)) '.']);