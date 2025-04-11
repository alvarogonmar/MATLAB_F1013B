%% create a matrix
close all;
clear;
clc;

% Define grid using linspace
x_vals = linspace(-5, 5, 100); % 100 point from -5 to 5
y_vals = linspace(-5, 5, 100);

% Preallocate the scalar field matrix
f = zeros(length(y_vals), length(x_vals));

% fill the matrix using nested for loops
for i = 1:length(y_vals)
    for j = 1:length(x_vals)
        x = x_vals(j);
        y = y_vals(i);
        f(i, j) = sin(x) * cos(y); % example scalar field
    end
end

disp(['Size of f: ', num2str(size(f,1)), ' rows x', num2str(size(f,2)), ' columns'])

% Create meshgrid for plotting
[X, Y] = meshgrid(x_vals, y_vals);

% Surface plot
figure;
surf(X, Y, f, 'EdgeColor', 'none');
colormap parula;
colorbar;
title('Surface Plot of Scalar Field');
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
hold on;

% view settings
view(45, 30);
grid on;
%% Part 2
% Define the grid
[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
