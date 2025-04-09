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

% Part 1: Compute the divergence and curl of the 2D vector field
% The 2D vector field is F(x, y) = (-y, x)
syms x y
F = [-y, x]; % Vector field

% Compute divergence: div(F) = ∂F1/∂x + ∂F2/∂y
div_F = divergence(F, [x, y]);
fprintf('Divergence of the field F(x, y): %s\n', div_F);

% Compute curl: curl(F) = ∂F2/∂x - ∂F1/∂y (in 2D)
curl_F = curl([F, 0], [x, y, 0]); % Extend to 3D for curl computation
fprintf('Curl of the field F(x, y): %s\n', curl_F(3));

% Part 2: Use meshgrid and quiver to plot the vector field
% Define the grid for x and y
[x_grid, y_grid] = meshgrid(-3:0.25:3, -3:0.25:3);

% Evaluate the components of the vector field
Fx = -y_grid; % F1 component
Fy = x_grid;  % F2 component

% Plot the vector field using quiver
figure;
quiver(x_grid, y_grid, Fx, Fy, 'b');
xlabel('x');
ylabel('y');
title('2D Vector Field F(x, y) = (-y, x)');
axis equal;
grid on;

% Part 3: Overlay streamlines using plot
% Define starting points for streamlines
startX = -3:0.5:3; % x-coordinates of starting points
startY = -3 * ones(size(startX)); % y-coordinates of starting points

% Plot the streamlines
hold on;
streamline(x_grid, y_grid, Fx, Fy, startX, startY);
legend('Vector Field', 'Streamlines');

% Part 4: Extend the problem to 3D
% Define a 3D vector field F(x, y, z)
syms z
F3D = [-y, x, z]; % Vector field in 3D
fprintf('3D Vector Field: F(x, y, z) = %s\n', F3D);

% Define the grid for x, y, and z
[x3D, y3D, z3D] = meshgrid(-3:1:3, -3:1:3, -3:1:3);

% Evaluate the components of the 3D vector field
Fx3D = -y3D;
Fy3D = x3D;
Fz3D = z3D;

% 3D quiver plot of the vector field
figure;
quiver3(x3D, y3D, z3D, Fx3D, Fy3D, Fz3D, 'r');
xlabel('x');
ylabel('y');
zlabel('z');
title('3D Vector Field F(x, y, z) = (-y, x, z)');
grid on;

% Part 5: Draw divergence of the 2D vector field
% Compute the divergence numerically
div_F_numeric = Fx ./ x_grid + Fy ./ y_grid;

% Plot the divergence as a scalar field
figure;
contourf(x_grid, y_grid, div_F_numeric, 20, 'LineColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('Divergence of the 2D Vector Field');
grid on;