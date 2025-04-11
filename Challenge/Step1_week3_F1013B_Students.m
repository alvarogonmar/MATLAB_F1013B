%%%%%%%%%%%%%%%%%%%%%%
%Step1_Session3_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%TO-DO%
%Using the linspace command, and the new template code included here,
%construct a "continuous" electric charge distribution, using a total
%electric charge Q per electrode, and a number of electric charges Nq to
%discretize it.

clear;      % Limpia todas las variables del espacio de trabajo.
clc;        % Limpia la ventana de comandos.
clf;        % Limpia todas las figuras abiertas.

%----------------- PARÁMETROS INICIALES ------------------%
Lp = 3.5;                      % Longitud del electrodo positivo en mm.
Ln = 2.5;                      % Longitud del electrodo negativo en mm.
t = 0.02;                      % Espesor adicional de los electrodos en mm.
d = 0.4;                       % Distancia entre los electrodos en mm.
p = 0.01;                      % Paso de discretización en mm.

% Características de los electrodos
ke = 1/(4*pi*8.85*10^-12);     % Constante de Coulomb, utilizada para calcular campos eléctricos.
Q = 1e-3;                      % Carga total por electrodo en Coulombs.
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

%--------------------- DISTRIBUCIÓN DE CARGAS ------------------%
% Cálculo del diferencial de carga
dq = Q / Nq;                   % Magnitud de cada carga discreta.

% Posiciones de las cargas positivas
yp = linspace(-(1-p)*Lp/2, (1-p)*Lp/2, Nq); % Posiciones Y de las cargas positivas.
xp(1:Nq) = -d/2 - t/2;                      % Posiciones X de las cargas positivas (constante).

% Posiciones de las cargas negativas
yn = linspace(-(1-p)*Ln/2, (1-p)*Ln/2, Nq); % Posiciones Y de las cargas negativas.
xn(1:Nq) = d/2 + t/2;                       % Posiciones X de las cargas negativas (constante).

%--------------------- VISUALIZACIÓN ------------------%
% Verificación de posicionamiento de las cargas
plot(xp, yp, '*', 'Color', colorP);         % Dibuja las posiciones de las cargas positivas.
hold on
plot(xn, yn, '*', 'Color', colorN);         % Dibuja las posiciones de las cargas negativas.

% Configuración de la gráfica
hold on
axis([xmin xmax ymin ymax])    % Ajusta los límites de los ejes.
xlabel('x position, mm')       % Etiqueta del eje x.
ylabel('y position, mm')       % Etiqueta del eje y.
title('Dielectrophoresis (No gradient)') % Título del gráfico.
grid on                        % Activa la cuadrícula.

% Dibujo de los electrodos
patch('Faces', facesP, 'Vertices', vertices2d, 'FaceColor', colorP); % Dibuja el electrodo positivo.
patch('Faces', facesN, 'Vertices', vertices2d, 'FaceColor', colorN); % Dibuja el electrodo negativo.