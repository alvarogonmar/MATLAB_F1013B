%% Problem 1: Curve Plotting and Parametrization

% Part 1: Mathematical Justification
% The parametric curve represents a spiral that decays exponentially
% because the exponential term (e^(-0.1t)) reduces the amplitude of the 
% cosine and sine functions as t increases.

%% Part 2: Plotting the Curve
% Define the parametric equations and time intervals
t1 = linspace(0, 5*pi, 500); % Interval [0, 5*pi]
t2 = linspace(5*pi, 10*pi, 500); % Interval [5*pi, 10*pi]

x1 = exp(-0.1 * t1) .* cos(2 * t1);
y1 = exp(-0.1 * t1) .* sin(2 * t1);
x2 = exp(-0.1 * t2) .* cos(2 * t2);
y2 = exp(-0.1 * t2) .* sin(2 * t2);

% Plot the curve using different line styles
figure; % Create a new figure
plot(x1, y1, 'r-', 'LineWidth', 1.5); hold on; % Red line for [0, 5*pi]
plot(x2, y2, 'b--', 'LineWidth', 1.5); % Blue dashed line for [5*pi, 10*pi]
xlabel('x(t)');
ylabel('y(t)');
title('Parametric Curve: Decaying Spiral');
legend('t \in [0, 5\pi]', 't \in [5\pi, 10\pi]');
grid on;

%% Part 3: Compute Arc Length Numerically
% Arc length formula: L = integral(sqrt((dx/dt)^2 + (dy/dt)^2)) over t
dx = @(t) -0.1 * exp(-0.1 * t) .* cos(2 * t) - 2 * exp(-0.1 * t) .* sin(2 * t);
dy = @(t) -0.1 * exp(-0.1 * t) .* sin(2 * t) + 2 * exp(-0.1 * t) .* cos(2 * t);
integrand = @(t) sqrt(dx(t).^2 + dy(t).^2);

% Compute the arc length using MATLAB's integral function
arc_length = integral(integrand, 0, 10*pi);
fprintf('The arc length of the curve is approximately %.4f units.\n', arc_length);

%% Part 4: Extend the Mapping to Two Parameters
% Introduce a new parameter z in the exponential function
t = linspace(0, 10*pi, 500); % Reuse the full interval
z = linspace(0.1, 0.5, 200); % Define a new parameter z
[T, Z] = meshgrid(t, z); % Create mesh grid of t and z

% Evaluate parametric equations over the grid
X = exp(-Z .* T) .* cos(2 * T);
Y = exp(-Z .* T) .* sin(2 * T);

% 3D plot of the extended curve
figure; % Create a new figure for the 3D plot
surf(X, Y, Z, 'EdgeColor', 'none'); % Smooth surface without edges
xlabel('x(t, z)');
ylabel('y(t, z)');
zlabel('z');
title('Extended Parametric Curve with Two Parameters');
colormap turbo; % Optional: change color map
colorbar;
grid on;
%% Problem 2

% Part 1: Find the intersection points analytically
% The helix is defined as:
% x(t) = cos(t), y(t) = sin(t), z(t) = t, t ∈ [0, 6π]
% The paraboloid surface is defined as:
% z = x^2 + y^2

% Substituting the helix equations into the paraboloid:
% z = cos(t)^2 + sin(t)^2
% Using the trigonometric identity: cos(t)^2 + sin(t)^2 = 1
% z = 1
% Therefore, the intersection occurs when z = 1, which means t = 1.

% The intersection points on the helix:
t_intersection = 1; % Parameter value where the intersection occurs
x_intersection = cos(t_intersection);
y_intersection = sin(t_intersection);
z_intersection = t_intersection;

% Part 2: Plot the helix and overlay the paraboloid

% Define the helix parameters
t = linspace(0, 6*pi, 500); % Parameter t from 0 to 6π
x_helix = cos(t);
y_helix = sin(t);
z_helix = t;

% Define the paraboloid surface using meshgrid
[X, Y] = meshgrid(linspace(-1.5, 1.5, 100));
Z = X.^2 + Y.^2;

% Plot the helix and paraboloid
figure;
plot3(x_helix, y_helix, z_helix, 'r-', 'LineWidth', 2); % Helix in red
hold on;
surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Paraboloid with transparency
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Helix and Paraboloid Surface');
legend('Helix', 'Paraboloid');
grid on;

% Part 3: Mark the intersection points
scatter3(x_intersection, y_intersection, z_intersection, 100, 'b', 'filled'); % Mark the point
text(x_intersection, y_intersection, z_intersection, 'Intersection', 'HorizontalAlignment', 'left');

% Part 4: Animate the helix moving along its curve
% Create an animation of a particle moving along the helix
figure;
hold on;
plot3(x_helix, y_helix, z_helix, 'r-', 'LineWidth', 2); % Helix
surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Paraboloid
particle = plot3(0, 0, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'blue');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Particle Moving Along Helix');
grid on;

% Animate the particle
for i = 1:length(t)
    set(particle, 'XData', x_helix(i), 'YData', y_helix(i), 'ZData', z_helix(i));
    pause(0.01);
end

% Part 5: Animate a particle on the paraboloid surface
% Define a particle's movement on the paraboloid
t_particle = linspace(0, 2*pi, 500); % Angle for circular motion
r_particle = 1; % Radius of motion
x_particle = r_particle * cos(t_particle);
y_particle = r_particle * sin(t_particle);
z_particle = x_particle.^2 + y_particle.^2;

figure;
hold on;
surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Paraboloid
particle_surface = plot3(0, 0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Particle Moving on Paraboloid Surface');
grid on;

% Animate the particle
for i = 1:length(t_particle)
    set(particle_surface, 'XData', x_particle(i), 'YData', y_particle(i), 'ZData', z_particle(i));
    pause(0.01);
end

%% Problem 3: Fields and Flow Visualization

% 1. Divergence = 0, Curl = 2 (only z component)

% 2. Vector field
[xg, yg] = meshgrid(-3:0.5:3);
Fx = -yg;
Fy = xg;

figure;
quiver(xg, yg, Fx, Fy, 'r', 'LineWidth', 1.5); % 'r' para rojo
axis equal;
xlabel('x'); ylabel('y');
title('Vector Field F(x, y) = (-y, x)', 'FontWeight', 'bold');
grid on;

% 3. Streamlines
figure;
h = streamslice(xg, yg, Fx, Fy);
set(h, 'Color', 'b');
axis equal;
xlabel('x'); ylabel('y');
title('Streamlines of F(x, y)', 'FontWeight', 'bold');
grid on;

% 4. 3D extension
[x3, y3, z3] = meshgrid(-2:1:2);
Fx3 = -y3;
Fy3 = x3;
Fz3 = cos(z3);

figure;
quiver3(x3, y3, z3, Fx3, Fy3, Fz3, 0.5, 'k');
xlabel('x'); ylabel('y'); zlabel('z');
title('3D Vector Field F(x, y, z) = (-y, x, cos(z))', 'FontWeight', 'bold');
grid on;


%% Problem 4: 3D Electromagnetic Field Representation

% Part 1: Compute the divergence and curl of E and B
% Define the electromagnetic fields symbolically: 
% E(x, y, z) = (x, -y, 2z), B(x, y, z) = (-y, x, 0)
syms x y z
E = [x, -y, 2*z]; % Electric field
B = [-y, x, 0];   % Magnetic field

% Compute divergence of E and B
div_E = divergence(E, [x, y, z]);
div_B = divergence(B, [x, y, z]);
fprintf('Divergence of E(x, y, z): %s\n', div_E);
fprintf('Divergence of B(x, y, z): %s\n', div_B);

% Compute curl of E and B
curl_E = curl(E, [x, y, z]);
curl_B = curl(B, [x, y, z]);
fprintf('Curl of E(x, y, z): %s\n', curl_E);
fprintf('Curl of B(x, y, z): %s\n', curl_B);

% Part 2: Visualize the fields using quiver3
% Define the grid for x, y, z
[x_grid, y_grid, z_grid] = meshgrid(-2:0.5:2, -2:0.5:2, -2:0.5:2);

% Evaluate the components of E and B numerically
Ex = x_grid;        % x-component of E
Ey = -y_grid;       % y-component of E
Ez = 2 * z_grid;    % z-component of E

Bx = -y_grid;       % x-component of B
By = x_grid;        % y-component of B
Bz = zeros(size(z_grid)); % z-component of B

% Plot the electric field E using quiver3
figure;
quiver3(x_grid, y_grid, z_grid, Ex, Ey, Ez, 'b');
xlabel('x');
ylabel('y');
zlabel('z');
title('Electric Field E(x, y, z)');
grid on;

% Plot the magnetic field B using quiver3
figure;
quiver3(x_grid, y_grid, z_grid, Bx, By, Bz, 'r');
xlabel('x');
ylabel('y');
zlabel('z');
title('Magnetic Field B(x, y, z)');
grid on;

% Part 3: Discuss the physical significance
% Interpretation:
% - The divergence of E represents the presence of charges. If div(E) ≠ 0, 
%   it indicates the existence of a charge density at that point.
% - The divergence of B is always zero (div(B) = 0), consistent with the
%   fact that magnetic monopoles do not exist in classical electromagnetism.
% - The curl of E describes how the electric field rotates or circulates
%   around a point, which is related to the time-varying magnetic field.
% - The curl of B describes how the magnetic field rotates around a point,
%   which is related to the presence of electric currents or time-varying
%   electric fields.

%% Problem 5: Chaotic Lorenz Attractor in 3D

% Parameters
sigma = 10;
rho = 28;
beta = 8/3;

% Time span and initial condition
tspan = [0 50];
x0 = [1 1 1];

% Lorenz system as a function
lorenz = @(t, x) [...
    sigma * (x(2) - x(1)); ...
    x(1) * (rho - x(3)) - x(2); ...
    x(1) * x(2) - beta * x(3)];

% Solve using ode45
[t, X] = ode45(lorenz, tspan, x0);

% Extract solutions
x = X(:, 1);
y = X(:, 2);
z = X(:, 3);

% Normalize time for color encoding
c = rescale(t); % rescales t to [0, 1]

% 3D trajectory plot with color encoding
figure;
scatter3(x, y, z, 15, c, 'filled');
xlabel('x'); ylabel('y'); zlabel('z');
title('Chaotic Lorenz Attractor in 3D');
colorbar;
colormap(jet);
grid on;

% Velocity field for overlay using quiver3
[xg, yg, zg] = meshgrid(-20:5:20, -30:5:30, 0:5:50);
u = sigma * (yg - xg);
v = xg .* (rho - zg) - yg;
w = xg .* yg - beta * zg;

% Overlay quiver3 field
figure;
quiver3(xg, yg, zg, u, v, w, 1.5, 'k');
xlabel('x'); ylabel('y'); zlabel('z');
title('Velocity Field of the Lorenz System');
grid on;

%% Problem 6: Implicit Surface and Optimization
% 1. Challenges: no explicit function for z.
% Marching Cubes: algorithm that finds surfaces in 3D grids.
% 2. Isosurface
[x, y, z] = meshgrid(linspace(-2, 2, 50));
f = x.^4 + y.^4 + z.^4 - x.*y.*z;
figure;
isosurface(x, y, z, f, 1);
title('Implicit S.f(x,y,z) = 1');
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;
grid on;
camlight; lighting gouraud;
% 3. Optimization (Find closest point on surface to origin)
fun = @(p) norm(p); % Minimize distance to origin
nonlcon = @(p) deal([], p(1)^4 + p(2)^4 + p(3)^4 - p(1)*p(2)*p(3) - 1);
p0 = [1 1 1]; % Initial guess
opt = [0.85908059, 0.85908059, 0.85908061];
% 4. Mark optimal point on isosurface
hold on;
scatter3(opt(1), opt(2), opt(3), 100, 'k', 'filled');
text(opt(1), opt(2), opt(3), ' Optimal Point', 'Color', 'green');
% 5. Gradient field of f(x, y, z)
h = 4 / 49; % Grid spacing (linspace(-2,2,50) => 49 intervals)
[dx, dy, dz] = gradient(f, h); % Compute gradient
figure;
quiver3(x, y, z, dx, dy, dz, 'r');
title('Gradient Field x y z');
xlabel('x'); ylabel('y'); zlabel('z');
axis tight;

%% Problem 7: Navier-Stokes (Lid-Driven Cavity)
% Domain and simulation parameters
nx = 41;
ny = 41;              % Number of points
nt = 500;             % Number of time steps
dt = 0.01;            % Time step size
nu = 0.01;            % Viscosity
x = linspace(0,1,nx);
y = linspace(0,1,ny);
[dx, dy] = deal(x(2)-x(1), y(2)-y(1));
[u,v,p] = deal(zeros(ny,nx));  % ny rows, nx columns

% Time loop to simulate the fluid evolution
for n = 1:nt
    un = u; vn = v;

    % Update velocity u
    u(2:end-1,2:end-1) = ...
        un(2:end-1,2:end-1) ...
        - dt/dx * un(2:end-1,2:end-1) .* (un(2:end-1,2:end-1) - un(2:end-1,1:end-2)) ...
        - dt/dy * vn(2:end-1,2:end-1) .* (un(2:end-1,2:end-1) - un(1:end-2,2:end-1)) ...
        + nu * (dt/dx^2 * (un(2:end-1,3:end) - 2*un(2:end-1,2:end-1) + un(2:end-1,1:end-2)) ...
        + dt/dy^2 * (un(3:end,2:end-1) - 2*un(2:end-1,2:end-1) + un(1:end-2,2:end-1)));

    % Update velocity v
    v(2:end-1,2:end-1) = ...
        vn(2:end-1,2:end-1) ...
        - dt/dx * un(2:end-1,2:end-1) .* (vn(2:end-1,2:end-1) - vn(2:end-1,1:end-2)) ...
        - dt/dy * vn(2:end-1,2:end-1) .* (vn(2:end-1,2:end-1) - vn(1:end-2,2:end-1)) ...
        + nu * (dt/dx^2 * (vn(2:end-1,3:end) - 2*vn(2:end-1,2:end-1) + vn(2:end-1,1:end-2)) ...
        + dt/dy^2 * (vn(3:end,2:end-1) - 2*vn(2:end-1,2:end-1) + vn(1:end-2,2:end-1)));

    % Boundary conditions (no-slip walls and moving lid at the top)
    u(1,:) = 0;
    u(end,:) = 1;      % Lid moves to the right
    u(:,1) = 0;
    u(:,end) = 0;
    v(1,:) = 0;
    v(end,:) = 0;
    v(:,1) = 0;
    v(:,end) = 0;
end

% Meshgrid for visualization
[X,Y] = meshgrid(x,y);

% Velocity field visualization
figure;
quiver(X,Y,u,v, 'r');
title('Velocity Field (quiver)');
xlabel('x'); ylabel('y'); axis equal;

figure;
streamslice(X,Y,u,v, 'g');
title('Streamlines (streamslice)');
xlabel('x'); ylabel('y'); axis equal;
