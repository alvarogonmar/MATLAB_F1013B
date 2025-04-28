%%%%%%%%%%%%%%%%%%%%%%
%Step0_Session4_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%Continue with the latest version of your code. Do the following.

%TO-DO%
%(1) Using now the concept of electric potential V(i,j) and its relation with the
%E-Field, calculate (using the same three nested for loops preferrably),
%the E field using Matlab's "gradient" function.  

%(2) Using again the command streamslice, create a new canvas to plot the
%aforementioned electric field, using the gradient.

%(3) Add and plot the equipotential lines using the "contour" (use Matlab documentation) command in
%both associated canvases.

%(4) Test the robustness of you code via changing the Nq, Lp and Ln
%values (what happens with the field lines?) - Expect questions of this
%sort in your Oral Exam.

clear;
clc;
clf;

%-----------------?------------------%
Lp=3.5;                      % Longitud placa positiva
Ln=2.5;                      % Longitud placa negativa
t=0.02;                      % Grosor de los electrones
d=0.4;                       % Distancia entre placas
p=0.01;                      % Distancia minima entre las cargas

% Define the characteristics of the electrodes
ke=1/(4*pi*8.85*10^-12);     % Constante del campo electrico
Q=1e-3;                      % Carga electrica en Coulombs
Nq=28;                       % Numero de cargas

%-----------------?------------------%

xmin_original=-d/2-3*t;  xmax_original=-xmin_original;  % 
xmin=-d/2-3*t;  xmax=-xmin;                             % 
ymin=2*(-Lp/2);   ymax=-ymin;                           % 

%-----------------Asegurar que las dimenciones que el area de trabajo sea el adecuado------------------%
if ymin <= -1  
    if xmin >= -0.5 && xmax <= 0.5
            xmin = -1.5;
            xmax = -xmin;
    end
end

%---------------% Se define una malla de puntos en el espacio 2D --------------------%
Ny=30;  Nx=Ny;
x=linspace(xmin, xmax, Nx); y=linspace(ymin, ymax, Ny); % Dibujar puntos en x y y

%-----------------Parches------------------%

vertices2d=[[-d/2-t,Lp/2]    %1 Vertice 1
    [-d/2,Lp/2]              %2 Vertice 2
    [-d/2,-Lp/2]             %3  Vertice 3
    [-d/2-t,-Lp/2]           %4 Vertice 4
    [d/2,Ln/2]               %5 Vertice 5
    [d/2+t,Ln/2]             %6  Vertice 6
    [d/2+t,-Ln/2]            %7  Vertice 7
    [d/2,-Ln/2]];            %8  Vertice 8

%  Las caras del cubo para conectar los vertices
facesP=[1 2 3 4 1];
facesN=[5 6 7 8 5];

% 
colorP=[0.95,0,0];           % Ponerle colores a las caras color rojo
colorN=[0,0,0.7];            % Color azul

%---------------------Start positioning charges------------------%
% Define a linear charge differential
dq=Q/Nq;                     % Charge differential magnitude

% Define the positions of the charges
yp=linspace(-(1-p)*Lp/2,(1-p)*Lp/2,Nq); % Positive charges Y positions, colorcar las cargas sobre las placas
xp(1:Nq)=-d/2-t/2;                      % Positive charges X positions
yn=linspace(-(1-p)*Ln/2,(1-p)*Ln/2,Nq);  % Negative charges Y positions
xn(1:Nq)= d/2+t/2;                                % Negative charges X positions

%Uncomment and complete the following template code to find out if
%you are placing you charges correctly. Additionally, play with some of the
%parameters above, to make sure that the dicrete charges are distributed
%between your plates. Finally, change the lenght of your plates, and see if
%your code responds accordingly.

plot(xp,yp,'*')
hold on
plot(xn,yn,'*') % Same but for negative charges


%-------Electric field calculation for every XY point (No gradient)-------%
% 
% % Initialize potential here
  V(1:Nx,1:Ny)=0;
% 
% % % Initialize electric field components
  % Inicializar componentes del campo eléctrico
Ex = zeros(Nx, Ny);  % Componente X del campo eléctrico
Ey = zeros(Nx, Ny);  % Componente Y del campo eléctrico
 
% % % Calculate electric field components
% %Three nested for loops start here...
 
  for i = 1:Nx % Loop for x coordinate
      for j = 1:Ny % Loop for y coordinate
          for k = 1:Nq % Loop for every charge

              rx = x(i);
              ry = y(j);
              rxp = xp(k); % Posiciones de x de las cargas positivas
              ryp = yp(k);
              rxn = xn(k);
              ryn = yn(k);

               % Electric Field
              r3psep = (sqrt((rx-rxp)^2 + (ry-ryp)^2))^3; % Distancia a la carga positiva
              r3nsep = (sqrt((rx - rxn)^2 + (ry - ryn)^2))^3; % Distancia a la carga negativa

              Ex(i, j) = Ex(i, j) + ke * dq * (rx - rxp) / r3psep - ke * dq * (rx - rxn)/r3nsep;% Componente x de E
              Ey(i, j) = Ey(i, j) + ke * dq * (ry - ryp) / r3psep - ke * dq * (ry -ryn)/r3nsep; % Componente y de E
               
                
              % Electric Potencial
              rpsep = (sqrt((rx-rxp)^2 + (ry-ryp)^2)); % Distancia a la carga positiva
              rnsep = (sqrt((rx - rxn)^2 + (ry - ryn)^2));
              V(i, j) = V(i, j) + ke * dq * (1/rpsep - 1/rnsep);
          end
      end
  end
 

% %Some place here you must calculate the E-field components.
% %Some place you must calculate V(i,j) potential for ith, jth position.

 
% %Three anidated for loops ends here...
% 
% % Calculate the electric field components using the gradient of potential
% % here
% 
[ExP, EyP] = gradient(V'); % Why do we need V' instead of V?
% Se usa la traspuesta V' porque en MATLAB las filas de una matriz representan el eje (y) 
% y las columnas el eje (x), lo cual es opuesto a la convención física donde columnas 
% corresponden a (x) y filas a (y).
ExP = -ExP;
EyP = -EyP;
% 
% % 
hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectrophoresis (No gradient)'
grid on
% 
% %Set aesthetics (for E-field plot) using the potential values recently added (just uncomment).
pcolor(x,y,V')                % Color map of the Voltage
colormap bone                 % Color
% 
% 
% %Use here the streamslice command to plot you calculated E-field...
streamslice(x,y,Ex',Ey',2);  % Electric field lines without gradient 
%                         % Why do we need Ex' and Ey' instead of Ex and Ey?
% 
% %Uncomment and comment (what are these new lines for?)
shading interp;
colorbar;
% 

patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);
 
% %%%-- Plot here now the E-Field using the gradient --%%%
% %Use previous lines as a template, use same aesthetics. 


figure;
hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectrophoresis (With gradient)'
grid on
% 
% %Set aesthetics (for E-field plot) using the potential values recently added (just uncomment).
pcolor(x,y,V')                % Color map of the Voltage
colormap bone

streamslice(x,y,ExP,EyP,2);
shading interp;
colorbar;
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);


%% Step 1

%-------------------------Erythrocytes configuration----------------------%
h = 0.12; % Define the characteristics of the erythrocytes

%                                 % Define an adjustable parameter (biology enters here) using the letter "h". Consider that 
                                   % h<0.03 will be associated with healthy blood and h>0.2 with infected blood 
                                   % (for a final set of given parameters).

rad=(xmax - xmin)/50;  % Define the erythrocyte radius as a fraction of the x space scale.

Ne = 10;                             % Create a variable to define the No. of erythrocytes to be modelled.

dt = 0.2;                             % Create a time step variable "dt" in Arbitrary units (use a value of 0.2)

qe = 1e-6;                             % Create a variable "qe" to define the erythrocyte electric charge 
                               % (equally positives and negatives - why?),
                               % use a value of 1 micro Coulombs.
                         
                            


%-------------------------Erythrocytes displacement-----------------------%
% Initialize the count of healthy and infected erythrocytes
healthy=0; infected=0;

% (3) Open a for loop that will run over the total number of erythrocytes (Ne) defined by the user, 
% and within it initialize the erythrocyte position (xe, ye) and velocities (Vx , Vy)... Think about it! :)

for ery = 1:Ne
    % (Step2 - Parte 1) Define a variable called path using animatedline
    path = animatedline('LineWidth',2,'LineStyle',':', 'Color','y'); % Traza la trayectoria del eritrocito
    
    % (Step2 - Parte 1) Introduce randomness to the erythrocyte's motion
    dx = h * rand() * rad;  % Desplazamiento aleatorio basado en el radio y parámetro h
    % Esto simula que cada eritrocito experimenta un pequeño "salto aleatorio" lateral
    fprintf('Erythrocyte: %d has an associated displacement of: %f\n', ery, dx);

    % Iniciar las posiciones x y y y velocidades
    xe = 0;
    ye = ymax;
    Vx = 0;
    Vy = -1;
    
    % (Step2 - Parte 2) Open while loop: mientras esté por encima de ymin
    while ye > ymin
        % Dibujar el movimiento del eritrocito en tiempo real
        addpoints(path, xe, ye); % Agrega el punto actual (xe, ye) a la trayectoria
        head = scatter(xe, ye, 100, 'filled', 'o', 'red'); % Representa la posición actual del eritrocito
        drawnow; % Actualiza el dibujo inmediatamente
        
        %--------------Force calculation-------------------%
        % (Step2 - Parte 2) Inicializar las componentes de la fuerza eléctrica
        Fx = 0;  % Componente X de la fuerza
        Fy = 0;  % Componente Y de la fuerza

        % (Dentro del while, abre un for para recorrer todas las cargas)
        for k = 1:Nq
            % Calcula las fuerzas de Coulomb desde cada carga positiva
            rxp = xp(k);  % X de la carga positiva
            ryp = yp(k);  % Y de la carga positiva
            r = sqrt((xe - rxp)^2 + (ye - ryp)^2); % Distancia a la carga positiva
            Fx = Fx + qe * ke * (xe - rxp) / r^3;  % Fuerza en X positiva
            Fy = Fy + qe * ke * (ye - ryp) / r^3;  % Fuerza en Y positiva

            % Calcula las fuerzas de Coulomb desde cada carga negativa
            rxn = xn(k);  % X de la carga negativa
            ryn = yn(k);  % Y de la carga negativa
            r = sqrt((xe - rxn)^2 + (ye - ryn)^2); % Distancia a la carga negativa
            Fx = Fx - qe * ke * (xe - rxn) / r^3;  % Fuerza en X negativa
            Fy = Fy - qe * ke * (ye - ryn) / r^3;  % Fuerza en Y negativa
        end

        % Actualizar velocidad y posición del eritrocito
        %----------------------Dynamics update----------------------%
        % Mass is assumed to be 1 kg
        ax = Fx;  % Acceleration in x
        ay = Fy;  % Acceleration in y
        
        Vx = Vx + Fx * dt;  % Velocidad en X
        Vy = Vy + Fy * dt;  % Velocidad en Y
        xe = xe + Vx * dt;  % Posición en X
        ye = ye + Vy * dt;  % Posición en Y
        
        % (Step2 - Parte 3) Después del for que recorre las cargas
        % Verificar si sigue arriba de ymin + dx
        if ye > ymin + dx
            delete(head); % Borra el punto anterior para no saturar la gráfica
        end
    end
end
