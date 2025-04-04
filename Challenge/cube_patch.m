clear;
clc;
close all;

% Define the vertices of the cube
vertices = [
    0 0 0; % v1
    1 0 0; % v2
    1 1 0; % v3
    0 1 0; % v4
    0 0 1; % v5
    1 0 1; % v6
    1 1 1; % v7
    0 1 1; % v8
];

% Define the faces of the cube using the vertex index
faces = [
  1 2 3 4; % Botton face
  5 6 7 8; % Top face
  1 2 6 5; % Front Face
  2 3 7 6; % Right face
  3 4 8 7; % Back Face
  4 1 5 8; % Left face
];

% Define colors
faceColors = [
    1 0 0; % Red
    0 1 0; % Green
    0 0 1; % Blue
    1 1 0; % Yellow
    1 0 1; % Magenta
    0 1 1 % Cyan
];

% Create figure and axes
figure;
ax = axes;
hold on;

% Draw cube faces with a single patch call
patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', faceColors, 'FaceColor', 'flat', 'EdgeColor', 'k');

% Set Axis properties
axis equal;
view(3);
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Cube using patch');