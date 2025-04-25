%%%%%%%%%%%%%%%%%%%%%%
%Step1_Session4_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%Now, you will work in configuring the Erythrocytes we are to simulate.

%TO-DO%
%(1) Go to bottom section of the code labelled "Erythrocytes configuration" and do
% as suggested there.     

%(2) Go to bottom section of the code labelled "Erythrocytes displacement" and do
% as suggested there.

%(3) Open a for loop that will run over the total number of erythrocytes (Ne) defined by the user, 
% and within it initialize the erythrocyte position (xe, ye) and velocities (Vx , Vy)... Think about it! :)
 
% (4) Next, inside this loop and after the initialization referred above, use (define) a while that will 
% run from the upmost Y-coordinate of the erythrocyte, say "ye", to the bottom-most (down-most)  
% Y-coordinate of the erythrocyte (say "ye>ymin"). 
% 
% (5)Finally, inside this while open another for loop that will run from charge 1 all the way to the last 
% charge in the plates (electrodes). 

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

%---------------------Start positioning charges------------------%
% Define a linear charge differential
%dq=?/?;                     % Charge differential magnitude

% Define the positions of the charges
yp=linspace(-(1-p)*Lp/2,(1-p)*Lp/2,Nq); % Positive charges Y positions
xp(1:Nq)=-d/2-t/2;                      % Positive charges X positions
yn=                                     % Negative charges Y positions
xn(?:?)=                                % Negative charges X positions

%Uncomment and complete the following template code to find out if
%you are placing you charges correctly. Additionally, play with some of the
%parameters above, to make sure that the dicrete charges are distributed
%between your plates. Finally, change the lenght of your plates, and see if
%your code responds accordingly.

%plot(<vector of positive charges X positions>, <vector of positive charges Y positions>,'*')
%hold on
%plot(?,?,'*') Same but for negative charges


%-------Electric field calculation for every XY point (No gradient)-------%

% % Initialize potential here
 V(?,?)=0;

% % Initialize electric field components
 Ex = zeros(?, ?);  %Where do you want this component to be calculated?
 Ey = zeros(?, ?);
 
% % Calculate electric field components
%Three nested for loops start here...

 for 




%Some place here you must calculate the E-field components.
%Some place you must calculate V(i,j) potential for ith, jth position.



 end

%Three anidated for loops ends here...

% Calculate the electric field components using the gradient of potential
% here

[?, ?] = gradient(V'); %Why do we need V' instead of V? 

% 
hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectrophoresis (No gradient)'
grid on

%Set aesthetics (for E-field plot) using the potential values recently added (just uncomment).
%pcolor(x,y,V')                % Color map of the Voltage
%colormap bone                 % Color


%Use here the streamslice command to plot you calculated E-field...
streamslice(?,?,?,?,?)  % Electric field lines without gradient 
                        % Why do we need Ex' and Ey' instead of Ex and Ey?

%Uncomment and comment (what are these new lines for?)
%shading interp;
%colorbar

% 
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);

%%%-- Plot here now the E-Field using the gradient --%%%
%Use previous lines as a template, use same aesthetics. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Step 1 Session4 Starts Here! %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------Erythrocytes configuration----------------------%
% Define the characteristics of the erythrocytes

%;                                 % Define an adjustable parameter (biology enters here) using the letter "h". Consider that 
                                   % h<0.03 will be associated with healthy blood and h>0.2 with infected blood 
                                   % (for a final set of given parameters).

%rad=(xmax_original-xmin_original)/40;  % Define the erythrocyte radius as a fraction of the x space scale.

%;                             % Create a variable to define the No. of erythrocytes to be modelled.

%;                             % Create a time step variable "dt" in Arbitrary units (use a value of 0.2)

%;                             % Create a variable "qe" to define the erythrocyte electric charge 
                               % (equally positives and negatives - why?),
                               % use a value of 1 micro Coulombs.
                         
                            


%-------------------------Erythrocytes displacement-----------------------%
% Initialize the count of healthy and infected erythrocytes
%healthy=0; infected=0;



% (3) Open a for loop that will run over the total number of erythrocytes (Ne) defined by the user, 
% and within it initialize the erythrocyte position (xe, ye) and velocities (Vx , Vy)... Think about it! :)
 
% (4) Next, inside this loop and after the initialization referred above, use (define) a while that will 
% run from the upmost Y-coordinate of the erythrocyte, say "ye", to the bottom-most (down-most)  
% Y-coordinate of the erythrocyte (say "ye>ymin"). 
% 
% (5)Finally, inside this while open another for loop that will run from charge 1 all the way to the last 
% charge in the plates (electrodes).


