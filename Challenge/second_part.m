clc; 
clear;
close all;

% Variables de configuración
v = 0.5; 
dx = -1; 
dx2 = 1;
dy = 0;  
dz = 0;  

largo = 3;  
largo2 = 3;  
ancho = -0.2;
alto = 1;

% Definir los vértices para el primer rectángulo (rojo)
vertices = [-v*ancho+dx  -v*largo+dy  -v+dz;  %v1
             -v*ancho+dx   v*largo+dy  -v+dz;  %v2
              v*ancho+dx   v*largo+dy  -v+dz;  %v3 
              v*ancho+dx  -v*largo+dy  -v+dz;  %v4
             -v*ancho+dx  -v*largo+dy   v+dz;  %v5
             -v*ancho+dx   v*largo+dy   v+dz;  %v6
              v*ancho+dx   v*largo+dy   v+dz;  %v7
              v*ancho+dx  -v*largo+dy   v+dz]; %v8
 
% Definir la cara del polígono usando los vértices definidos previamente
caras = [  
    1 2 3 4; % Cara inferior
    5 6 7 8; % Cara superior
    1 2 6 5; % Cara frontal
    3 4 8 7; % Cara posterior
    1 4 8 5; % Cara izquierda
    2 3 7 6]; % Cara derecha
 
% Definir los vértices para el segundo rectángulo (azul)
vertices2 = [-v*ancho+dx2  -v*largo+dy  -v+dz;  %v1
             -v*ancho+dx2   v*largo+dy  -v+dz;  %v2
              v*ancho+dx2   v*largo+dy  -v+dz;  %v3 
              v*ancho+dx2  -v*largo+dy  -v+dz;  %v4
             -v*ancho+dx2  -v*largo+dy   v+dz;  %v5
             -v*ancho+dx2   v*largo+dy   v+dz;  %v6
              v*ancho+dx2   v*largo+dy   v+dz;  %v7
              v*ancho+dx2  -v*largo+dy   v+dz]; %v8
         
% Definir la cara del polígono usando los vértices definidos previamente
caras2 = [  
    1 2 3 4; % Cara inferior
    5 6 7 8; % Cara superior
    1 2 6 5; % Cara frontal
    3 4 8 7; % Cara posterior
    1 4 8 5; % Cara izquierda
    2 3 7 6]; % Cara derecha

% Graficar las placas
hold on;
patch('Vertices', vertices, 'Faces', caras, 'FaceColor', 'r', 'EdgeColor', 'k'); % Placa roja
patch('Vertices', vertices2, 'Faces', caras2, 'FaceColor', 'b', 'EdgeColor', 'k'); % Placa azul

axis([-2 2 -2 2 -2 2]); % 
view(30,30); % 

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@ Principio de Superposición @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% Creación de la malla
x0 = -2;
x1 = 2;
y0 = -2;
y1 = 2;  
nx = 0.1; % Resolución más fina para mejorar la visualización
ny = 0.1;
       
[xGrid, yGrid] = meshgrid(x0:nx:x1, y0:ny:y1); % Parametrización de la malla
zGrid = xGrid * 0; % Creamos el parámetro z, pero lo mantenemos en cero

% Graficar la malla
plot(xGrid,yGrid,'.b');
h = quiver3(xGrid, yGrid, zGrid, xGrid, yGrid, zGrid, 'g');

hold on;

% Primer carga (Placa azul)
Qn1 = -20; % Carga negativa en la superficie azul
eps0 = 8.854e-12;  % Constante dieléctrica del medio
k = 1/(4*pi*eps0); % Constante eléctrica
Rxn1 = xGrid - (1); % Distancia en x desde la carga
Ryn1 = yGrid - (-largo2/2); % Distancia en y desde la carga
Rzn1 = zGrid; % Distancia en z desde la carga

% Cálculo del campo de una carga
R1 = sqrt(xGrid.^2 + yGrid.^2 + zGrid.^2);
Ex1 = k * Qn1 .* Rxn1 ./ R1.^3;
Ey1 = k * Qn1 .* Ryn1 ./ R1.^3;
Ez1 = k * Qn1 .* Rzn1 ./ R1.^3;

% Normalizar los componentes del vector dividiendo por E
i = Ex1 ./ E; 
j = Ey1 ./ E;
k = Ez1 ./ E; 

% Graficar
figure;
quiver3(xGrid, yGrid, zGrid, i, j, k, 'b'); % Flechas verdes para el campo

% Configuración final de la visualización
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Campo Eléctrico de la Primera Carga (Placa Azul)');
grid on;
hold off;

% Segunda carga (Placa roja)
Qp1 = 20; 
Rxp1 = xGrid - (-1); % Distancia en x desde la carga
Ryp1 = yGrid - (-1.5); % Distancia en y desde la carga
Rzp1 = zGrid; % Distancia en z desde la carga

% Cálculo del campo de una carga
R2 = sqrt(Rxp1.^2 + Ryp1.^2 + Rzp1.^2) + 1e-9; % Evitar división por cero

Ex2 = k * Qp1 .* Rxp1 ./ R2.^3;
Ey2 = k * Qp1 .* Ryp1 ./ R2.^3;
Ez2 = k * Qp1 .* Rzp1 ./ R2.^3;

% Aplicar el principio de superposición: Sumar los campos de las cargas anteriores
Ex_total = Ex1 + Ex2; 
Ey_total = Ey1 + Ey2;
Ez_total = Ez1 + Ez2; 

% Calcular la magnitud total del campo
E_total = sqrt(Ex_total.^2 + Ey_total.^2 + Ez_total.^2) + 1e-9; % Evitar división por cero

% Normalizar los componentes del vector dividiendo por E_total
i = Ex_total ./ E_total; 
j = Ey_total ./ E_total;
k = Ez_total ./ E_total; 

% Graficar
figure;
quiver3(xGrid, yGrid, zGrid, i, j, k, 'g'); % Flechas verdes para el campo

% Configuración de la visualización
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Campo Eléctrico de la Segunda Carga');
grid on;
axis equal;
view(30, 30);

% Tercera carga (Carga negativa 2)
Qn2 = -20; 
Rxn2 = xGrid - (1); 
Ryn2 = yGrid - (largo2 / 2); 
Rzn2 = zGrid - 0.5; 

R2 = sqrt(Rxn2.^2 + Ryn2.^2 + Rzn2.^2) + 1e-9; 

Ex2 = k * Qn2 .* Rxn2 ./ R2.^3;
Ey2 = k * Qn2 .* Ryn2 ./ R2.^3;
Ez2 = k * Qn2 .* Rzn2 ./ R2.^3;

Ex_total = Ex1 + Ex2; 
Ey_total = Ey1 + Ey2;
Ez_total = Ez1 + Ez2;

E_total = sqrt(Ex_total.^2 + Ey_total.^2 + Ez_total.^2) + 1e-9; 

i = Ex_total ./ E_total; 
j = Ey_total ./ E_total;
k = Ez_total ./ E_total; 

figure;
quiver3(xGrid, yGrid, zGrid, i, j, k, 'r'); 

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Campo Eléctrico de la Tercera Carga');
grid on;
axis equal;
view(30, 30);

% Cuarta carga (Carga positiva 2)
Qp2 = 20; 
Rxp2 = xGrid - (-1); 
Ryp2 = yGrid - (1.5); 
Rzp2 = zGrid - 0.5; 

R2 = sqrt(Rxp2.^2 + Ryp2.^2 + Rzp2.^2) + 1e-9;

Ex2 = k * Qp2 .* Rxp2 ./ R2.^3;
Ey2 = k * Qp2 .* Ryp2 ./ R2.^3;
Ez2 = k * Qp2 .* Rzp2 ./ R2.^3;

Ex_total = Ex1 + Ex2 + Ex3; 
Ey_total = Ey1 + Ey2 + Ey3;
Ez_total = Ez1 + Ez2 + Ez3;

E_total = sqrt(Ex_total.^2 + Ey_total.^2 + Ez_total.^2) + 1e-9; 

i = Ex_total ./ E_total; 
j = Ey_total ./ E_total;
k = Ez_total ./ E_total; 

figure;
quiver3(xGrid, yGrid, zGrid, i, j, k, 'b'); 

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Campo Eléctrico de la Cuarta Carga');
grid on;
axis equal;
view(30, 30);
