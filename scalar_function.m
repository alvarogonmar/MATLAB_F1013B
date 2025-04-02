% scalar function z = f(x,y)
clear;
clc;
close all;

% define the values (x, y) (mesh)
[x, y] = meshgrid(-2:0.05:2, -2:0.05:2);
% define the scalar function z = f(x,y)
z = x.* exp (-x.^2 - y.^2);
% plot the scalar functions along with the coordinates
surf(x, y, z);