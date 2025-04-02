clear;
clc;
close all;

% Define the domain
t = linspace(0, 2*pi, 1000);

% Define the mapping
x = cos(t);
y = sin(t);
z = (t);

figure;


for k = 1:10:length(t)
    % Clear figure
    clf;
    
    % Get data for instant k
    t_k = t(k);
    x_k = x(k);
    y_k = y(k);
    z_k = z(k);
    
    
    % Plot the current position of the particle
    plot3(x_k, y_k, z_k, '.r', 'MarkerSize', 100);
    hold on;
    grid on;

    % Set labels and title

    plot3(x,y,z, 'k')
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Particle at t = ' num2str(t_k) ' seconds']);   
    
    % Draw the current frame
    drawnow;

    % save the animation frames
    % vectorAnimation(end+1) = getframe(gcf);  
end

% Create video writer object
% miAnimacion = VideoWriter("circle_animation", "MPEG-4");
% miAnimacion.FrameRate = 20;

% Open the VideoWriter object
% open(miAnimacion);
% writeVideo(miAnimacion, vectorAnimation);
% close(miAnimacion);