% mesh grid
% first clear the variables of memory and workspace
clear;
clc;
close all;

m = 3; % define the size of the space
dm = 0.5; % define the increment of points

x = -m:dm:m; % x goes from -m to m
y = -m:dm:m; % in increments of dm

% get the mesh grid
[xx, yy] = meshgrid(x, y);

%plot the grid
plot(xx, yy);
grid on;