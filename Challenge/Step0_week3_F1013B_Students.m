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


% Limpieza de variables, ventana de comandos y figuras
clear;      % Limpia todas las variables del espacio de trabajo.
clc;        % Limpia la ventana de comandos.
clf;        % Limpia todas las figuras abiertas.

%----------------- PARÁMETROS INICIALES ------------------%
Lp = 3.5;                     % Longitud del electrodo positivo en mm.
Ln = 2.5;                     % Longitud del electrodo negativo en mm.
t = 0.02;                     % Espesor adicional de los electrodos en mm.
d = 0.4;                      % Distancia entre los electrodos en mm.
p = 0.01;                     % Paso de discretización en mm.

% Características de los electrodos
ke = 1/(4*pi*8.85*10^-12);     % Constante de Coulomb, utilizada para calcular campos eléctricos.
Q = 1e-3;                      % Magnitud de la carga en Coulombs.
Nq = 28;                       % Número de puntos discretos de carga.

%----------------- DEFINICIÓN DEL DOMINIO ------------------%
xmin_original = -d/2 - 3*t;    % Límite inferior original del eje x.
xmax_original = -xmin_original; % Límite superior original del eje x.
xmin = -d/2 - 3*t;             % Límite ajustado del eje x.
xmax = -xmin;                  % Límite ajustado del eje x.
ymin = 2*(-Lp/2);              % Límite inferior del eje y.
ymax = -ymin;                  % Límite superior del eje y.

% Bloque condicional para ajustar los límites si se cumplen ciertas condiciones
if ymin <= -1  
    if xmin >= -0.5 && xmax <= 0.5
        xmin = -1.5;           % Ajusta el límite inferior del eje x a -1.5.
        xmax = -xmin;          % Ajusta el límite superior del eje x a 1.5.
    end
end

% Discretización del dominio
Ny = 30;                       % Número de puntos en el eje y.
Nx = Ny;                       % Número de puntos en el eje x (iguales a Ny).
x = linspace(xmin, xmax, Nx);  % Discretización del eje x en Nx puntos.
y = linspace(ymin, ymax, Ny);  % Discretización del eje y en Ny puntos.

%----------------- GEOMETRÍA DE LOS ELECTRODOS ------------------%
vertices2d = [                 % Coordenadas de los vértices de los electrodos en 2D.
    [-d/2-t, Lp/2]             % Vértice 1 (electrodo positivo)
    [-d/2, Lp/2]               % Vértice 2
    [-d/2, -Lp/2]              % Vértice 3
    [-d/2-t, -Lp/2]            % Vértice 4
    [d/2, Ln/2]                % Vértice 5 (electrodo negativo)
    [d/2+t, Ln/2]              % Vértice 6
    [d/2+t, -Ln/2]             % Vértice 7
    [d/2, -Ln/2]               % Vértice 8
];

facesP = [1 2 3 4 1];          % Conexión de vértices para el electrodo positivo.
facesN = [5 6 7 8 5];          % Conexión de vértices para el electrodo negativo.
colorP = [0.95, 0, 0];         % Color del electrodo positivo (rojo).
colorN = [0, 0, 0.7];          % Color del electrodo negativo (azul).

%----------------- CONFIGURACIÓN DE LA GRÁFICA ------------------%
hold on                        % Permite añadir varios elementos a la gráfica.
axis([xmin xmax ymin ymax])    % Ajusta los límites de los ejes.
xlabel('x position, mm')       % Etiqueta del eje x.
ylabel('y position, mm')       % Etiqueta del eje y.
title('Dielectrophoresis (No gradient)') % Título del gráfico.
grid on                        % Activa la cuadrícula.

% Dibujo de los electrodos
patch('Faces', facesP, 'Vertices', vertices2d, 'FaceColor', colorP); % Dibuja el electrodo positivo.
patch('Faces', facesN, 'Vertices', vertices2d, 'FaceColor', colorN); % Dibuja el electrodo negativo.

