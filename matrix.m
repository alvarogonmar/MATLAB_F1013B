%% create a matrix
close all;
clear;
clc;

% Define grid using linspace
x_vals = linspace(-5, 5, 100); % 100 point from -5 to 5
y_vals = linspace(-5, 5, 100);

% Preallocate the scalar field matrix
f = zeros(length(y_vals), length(x_vals));

