% scalar function z = f(x,y)
clear;
clc;
close all;

% define the values (x, y) (mesh)
[x, y] = meshgrid(-3*pi:0.05:3*pi, -3*pi:0.05:3*pi);
% define the scalar function z = f(x,y)
z = sin(sqrt(x.^2 + y.^2));
% plot the scalar functions along with the coordinates
surf(x, y, z);