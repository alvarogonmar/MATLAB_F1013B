clear;
clc;
close all;

% define the values (x,y) (mesh)
[x, y] = meshgrid(-2:0.5:2, -2:0.5:2);

% define
fx = x;
fy = y;

% graficar con quiver
quiver(x, y, fx, fy, 'b')