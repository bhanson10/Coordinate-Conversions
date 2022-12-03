% Creating Inputs for Numerical Integration
Y0 = [20000; 0; 0; 0; 2.9; 1.8]; % [x; y; z; vx; vy; vz] [km, km/s]
tspan = [0 12*60*60]; % 12 hours [s]
options = odeset('RelTol', 1e-13); % Setting a tolerance
% Numerical Integration
[t, Y] = ode113(@ODE2BP, tspan, Y0, options);

% Pulling Position Data from Output
x = Y(:, 1); % [km]
y = Y(:, 2); % [km]
z = Y(:, 3); % [km]
xdot = Y(:, 4); % [km/s]
ydot = Y(:, 5); % [km/s]
zdot = Y(:, 6); % [km/s]
% % Setting up the Plot
figure; hold on; grid on;
ax = gca;
ax.FontSize = 16;
title('Velocity over Time - Cartesian', 'Interpreter', 'Latex', 'FontSize', 24)
xlabel('t (s)', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('Velocity (km/s)', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0, 43200])
ylim([-10 10])

% Create file name variable
filename = 'vel-cart.gif';

% Plotting the first iteration
cur_x_vel_plot = scatter(t(1), xdot(1), 'filled', 'b', 'DisplayName', 'X Velocity');
cur_y_vel_plot = scatter(t(1), ydot(1), 'filled', 'g', 'DisplayName', 'Y Velocity');
cur_z_vel_plot = scatter(t(1), zdot(1), 'filled', 'r', 'DisplayName', 'Z Velocity');

x_vel_plot = plot(t(1), x(1), 'b', 'LineWidth', 2, 'HandleVisibility','off');
y_vel_plot = plot(t(1), y(1), 'g', 'LineWidth', 2, 'HandleVisibility','off');
z_vel_plot = plot(t(1), z(1), 'r', 'LineWidth', 2, 'HandleVisibility','off');

lgd = legend;
lgd.Location = 'northeast';
pause(5);
% Iterating through the length of the time array
for k = 1:length(t)-1
    
    % Updating the point
    cur_x_vel_plot.XData = t(k);
    cur_x_vel_plot.YData = xdot(k);
    cur_y_vel_plot.XData = t(k);
    cur_y_vel_plot.YData = ydot(k);
    cur_z_vel_plot.XData = t(k);
    cur_z_vel_plot.YData = zdot(k);

    % Updating the lines
    x_vel_plot.XData = t(1:k);
    x_vel_plot.YData = xdot(1:k);
    y_vel_plot.XData = t(1:k);
    y_vel_plot.YData = ydot(1:k);
    z_vel_plot.XData = t(1:k);
    z_vel_plot.YData = zdot(1:k);
   
    % Delay
    pause(0.001)
    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime', (t(k+1)-t(k))/2500);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime', (t(k+1)-t(k))/2500);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dYdt = ODE2BP(t, Y)
    mu = 3.986*10^5; % Earth's gravitational parameter [km^3/s^2]
    x = Y(1); % [km]
    y = Y(2); % [km]
    z = Y(3); % [km]
    vx = Y(4); % [km/s]
    vy = Y(5); % [km/s]
    vz = Y(6); % [km/s]
    xddot = -mu/(x^2+y^2+z^2)^(3/2)*x; % [km/s^2]
    yddot = -mu/(x^2+y^2+z^2)^(3/2)*y; % [km/s^2]
    zddot = -mu/(x^2+y^2+z^2)^(3/2)*z; % [km/s^2]
    dYdt = [vx;vy;vz;xddot;yddot;zddot]; % Y'
end