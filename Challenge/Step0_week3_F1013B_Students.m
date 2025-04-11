%%%%%%%%%%%%%%%%%%%%%%
%Step0_Session3_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%This template code is very similar to the one you did in the previous 
%challenge session, with some additional lines of code that  
%will warranty that the visualized (plotted) E-field will look good, independently of us
%using the superposition principle to add all contributions to the total
%E-field at each point, or using the gradient of the electric potential to
%calculate the total E-field (also at each point).

%TO-DO%
%Investigate, run, and comment each line, identify the new block(s) of code and
%find out what they do (comment!), this will be part of your Deliverable No. 2.


clear;
clc;
clf;

%-----------------?------------------%
Lp=3.5;                      % 
Ln=2.5;                      % 
t=0.02;                      % 
d=0.4;                       % 
p=0.01;                      % 

% Define the characteristics of the electrodes
ke=1/(4*pi*8.85*10^-12);     % 
Q=1e-3;                      % What is this?
Nq=28;                       % And this?

%-----------------?------------------%

xmin_original=-d/2-3*t;  xmax_original=-xmin_original;  % 
xmin=-d/2-3*t;  xmax=-xmin;                             % 
ymin=2*(-Lp/2);   ymax=-ymin;                           % 

%-----------------?------------------%
if ymin <= -1  
    if xmin >= -0.5 && xmax <= 0.5
            xmin = -1.5;
            xmax = -xmin;
    end
end

%-----------------?------------------%
Ny=30;  Nx=Ny;
x=linspace(xmin, xmax, Nx); y=linspace(ymin, ymax, Ny);

%-----------------?------------------%

vertices2d=[[-d/2-t,Lp/2]    %1
    [-d/2,Lp/2]              %2
    [-d/2,-Lp/2]             %3  
    [-d/2-t,-Lp/2]           %4
    [d/2,Ln/2]               %5
    [d/2+t,Ln/2]             %6  
    [d/2+t,-Ln/2]            %7  
    [d/2,-Ln/2]];            %8  

%  
facesP=[1 2 3 4 1];
facesN=[5 6 7 8 5];

% 
colorP=[0.95,0,0];           % 
colorN=[0,0,0.7];            % 


% 

hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectrophoresis (No gradient)'
grid on


% 
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);


