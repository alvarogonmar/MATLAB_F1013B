%% How to calculate a derivative in Matlab

clear;
clc;
close all;

% define an anonymous function
f = @(x) x.^2 - 3*x - 4;

% First create a symbolic variable (since we are deriving with respect x) 
syms x;
