% Creating Inputs for Numerical Integration
a_d0 = 10^-5; mu = 3.986*10^5; a_d = [0 a_d0 0]; %a_d0 was varied based on the problem settings
Y0 = [20000; 0; 0; 0; 2.9; 1.8]; % [x; y; z; vx; vy; vz] [km, km/s]
tspan = [0 18*60*60]; % 12 hours [s]
options = odeset('RelTol', 1e-13); % Setting a tolerance
% Numerical Integration
[t, Y] = ode45(@(t, Y) ODE3BP(t, Y, mu, a_d), tspan, Y0, options);

% Pulling Position Data from Output
x = Y(:, 1); % [km]
y = Y(:, 2); % [km]
z = Y(:, 3); % [km]
xdot = Y(:, 4); % [km/s]
ydot = Y(:, 5); % [km/s]
zdot = Y(:, 6); % [km/s]

% % Setting up the Plot
figure; hold on
title('Two-Body Trajectory - Cartesian', 'Interpreter', 'Latex')
xlabel('x', 'Interpreter', 'Latex')
ylabel('y', 'Interpreter', 'Latex')
zlabel('z', 'Interpreter', 'Latex')
axis equal
grid minor
view(30, 30)

% Create file name variable
filename = '3bp-cart-2.gif';

% Creating/Plotting Spherical Earth
rm = 1000; % Radius of Earth [km]
[xEarth, yEarth, zEarth] = sphere(25);
surf(rm*xEarth,rm*yEarth,rm*zEarth, 'FaceColor', [0 1 0], 'HandleVisibility','off');

% Plotting the first iteration
p = plot3(x(1),y(1),z(1),'b', 'HandleVisibility','off');
m = scatter3(x(1),y(1),z(1),'filled','b', 'HandleVisibility','off');
pos_x = [0 x(1)]; pos_y = [0 y(1)]; pos_z = [0 z(1)];
vel_x = [x(1) (1000*xdot(1) + x(1))]; vel_y = [y(1) (1000*ydot(1) + y(1))]; vel_z = [z(1) (1000*zdot(1) + z(1))];
pos = plot3(pos_x, pos_y, pos_z, 'g', 'DisplayName', 'Position Vector');
vel = plot3(vel_x, vel_y, vel_z, 'r', 'DisplayName', 'Scaled Velocity Vector');
xlim([-20000 25000])
ylim([-1.5*10^4 1.5*10^4])
zlim([-1.5*10^4 1.5*10^4])
lgd = legend;
lgd.Location = 'northeast';

pause(5);
% Iterating through the length of the time array
for k = 1:length(t)-1
    
    % Updating the line
    p.XData = x(1:k);
    p.YData = y(1:k);
    p.ZData = z(1:k);
    % Updating the point
    m.XData = x(k); 
    m.YData = y(k);
    m.ZData = z(k);
    % Updating the Position vector
    pos_x = [0 x(k)]; pos_y = [0 y(k)]; pos_z = [0 z(k)];
    pos.XData = pos_x;
    pos.YData = pos_y;
    pos.ZData = pos_z;
    % Updating the Velocity vector
    vel_x = [x(k) (1000*xdot(k) + x(k))]; vel_y = [y(k) (1000*ydot(k) + y(k))]; vel_z = [z(k) (1000*zdot(k) + z(k))];
    vel.XData = vel_x;
    vel.YData = vel_y;
    vel.ZData = vel_z;
    % Updating the title
    title(sprintf('Two-Body Trajectory - Cartesian\nTime: %0.2f secs', t(k)),...
    'Interpreter','Latex');
    % Delay
    pause(0.001)
    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',(t(k+1)-t(k))/25000);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',(t(k+1)-t(k))/25000);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dYdt = ODE3BP(t, Y, mu, a_d)
    x = Y(1); % [km]
    y = Y(2); % [km]
    z = Y(3); % [km]
    vx = Y(4); % [km/s]
    vy = Y(5); % [km/s]
    vz = Y(6); % [km/s]
    xddot = -mu/(x^2+y^2+z^2)^(3/2)*x + a_d(1); % [km/s^2]
    yddot = -mu/(x^2+y^2+z^2)^(3/2)*y + a_d(2); % [km/s^2]
    zddot = -mu/(x^2+y^2+z^2)^(3/2)*z + a_d(3); % [km/s^2]
    dYdt = [vx;vy;vz;xddot;yddot;zddot]; % Y'
end