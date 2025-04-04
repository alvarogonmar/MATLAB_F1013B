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

