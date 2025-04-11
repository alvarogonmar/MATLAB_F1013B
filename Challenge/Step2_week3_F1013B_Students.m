%%%%%%%%%%%%%%%%%%%%%%
%Step2_Session3_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%TO-DO%
%(1) Using three nested for loops, one from i=1:Nx, one from j=1:Ny to cover space, and one from k=1:Nq to cover 
%all charges in the plates, calculate the total electric field inside and around the parallel plates.
%Use the hint code provided, and remember that a single point in space has contributions from each and all 
%electric charges in the plates (i.e. all positive and all negative charges).

%(2) Using a new Matlab command called "streamslice", plot the aforementioned electric field, you already 
% know how does it has to look.

%(3) Test the robustness of you code via changing the Nq, Lp and Ln
%values (what happens with the field lines?) - Expect questions of this
%behavior in your Oral Exam.


clear;
clc;
clf;

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

% Ajustar los límites del dominio, si es necesario
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

%----------------- DISTRIBUCIÓN DE CARGAS ------------------%
% Cálculo del diferencial de carga
dq = Q / Nq;                   % Magnitud de cada carga discreta.

% Posiciones de las cargas positivas
yp = linspace(-(1-p)*Lp/2, (1-p)*Lp/2, Nq); % Posiciones Y de las cargas positivas.
xp(1:Nq) = -d/2 - t/2;                      % Posiciones X de las cargas positivas (constante).

% Posiciones de las cargas negativas
yn = linspace(-(1-p)*Ln/2, (1-p)*Ln/2, Nq); % Posiciones Y de las cargas negativas.
xn(1:Nq) = d/2 + t/2;                       % Posiciones X de las cargas negativas (constante).

% Verificación gráfica de las posiciones de las cargas
plot(xp, yp, '*', 'Color', [0.95, 0, 0]);   % Dibuja las posiciones de las cargas positivas.
hold on;
plot(xn, yn, '*', 'Color', [0, 0, 0.7]);    % Dibuja las posiciones de las cargas negativas.

%----------------- CÁLCULO DEL CAMPO ELÉCTRICO ------------------%
% Inicializar componentes del campo eléctrico
Ex = zeros(Nx, Ny); % Componente x del campo
Ey = zeros(Nx, Ny); % Componente y del campo

% Calcular componentes del campo eléctrico mediante tres bucles anidados
for i = 1:Nx % Recorre las posiciones en el eje x
    for j = 1:Ny % Recorre las posiciones en el eje y
        for k = 1:Nq % Recorre cada carga discreta
            % Contribución de las cargas positivas
            rxp = x(i) - xp(k); % Diferencia en x
            ryp = y(j) - yp(k); % Diferencia en y
            rp = sqrt(rxp^2 + ryp^2)^3; % Distancia al cubo
            Ex(i, j) = Ex(i, j) + ke * (dq) * rxp / rp; % Componente x
            Ey(i, j) = Ey(i, j) + ke * (-dq) * ryp / rp; % Componente y

            % Contribución de las cargas negativas
            rxn = x(i) - xn(k); % Diferencia en x
            ryn = y(j) - yn(k); % Diferencia en y
            rn = sqrt(rxn^2 + ryn^2)^3; % Distancia al cubo
            Ex(i, j) = Ex(i, j) + ke * dq * rxn / rn; % Componente x
            Ey(i, j) = Ey(i, j) + ke * dq * ryn / rn; % Componente y
        end
    end
end

%----------------- VISUALIZACIÓN DEL CAMPO ELÉCTRICO ------------------%
hold on;
axis([xmin xmax ymin ymax]);
xlabel('x position, mm');
ylabel('y position, mm');
title('Dielectrophoresis (No gradient)');
grid on;

% Visualizar líneas de campo eléctrico usando "streamslice"
streamslice(x, y, Ex', Ey', 2); % Las líneas de campo eléctrico
                                % Nota: Se transponen Ex y Ey para ajustarse
                                % al formato esperado por streamslice.

% Dibujo de los electrodos
patch('Faces', facesP, 'Vertices', vertices2d, 'FaceColor', [0.95, 0, 0]); % Dibuja el electrodo positivo
patch('Faces', facesN, 'Vertices', vertices2d, 'FaceColor', [0, 0, 0.7]); % Dibuja el electrodo negativo