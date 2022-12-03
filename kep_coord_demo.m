% Creating Inputs for Numerical Integration
Y0 = [20000; 0; 0; 0; 2.9; 1.8]; % [x; y; z; vx; vy; vz] [km, km/s] (MEO)
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

mu = 3.986*10^5; % Earth's gravitational parameter [km^3/s^2]
[a, mag_e, e_vec, inc, O, w, f] = coordinate_conversion(x, y, z, xdot, ydot, zdot, mu);

% % Setting up the Plot
figure; hold on
title('Two-Body Trajectory - Keplerian', 'Interpreter', 'Latex')
xlabel('x', 'Interpreter', 'Latex')
ylabel('y', 'Interpreter', 'Latex')
zlabel('z', 'Interpreter', 'Latex')
axis equal
grid minor
xlim([-10000 26000])
ylim([-1.5*10^4 1.5*10^4])
zlim([-1.5*10^4 1.5*10^4])
view(30, 30)

% Create file name variable
filename = '2bp-kep.gif';

% Creating/Plotting Spherical Earth
rm = 1000; % Radius of Earth [km]
[xEarth, yEarth, zEarth] = sphere(25);
surf(rm*xEarth,rm*yEarth,rm*zEarth, 'FaceColor', [0 1 0], 'HandleVisibility','off');

% Plotting the first iteration

% Reference Plane
ref_x = [24000 24000 -10000 -10000];
ref_y = [-1.4*10^4 1.4*10^4 1.4*10^4 -1.4*10^4];
ref_z = [0 0 0 0];
patch(ref_x,ref_y,ref_z, 'k', 'DisplayName', 'Reference Plane')
hold off
alpha(0.4)

% Trajectory Plane
hold on
patch(x(1:160),y(1:160), z(1:160), 'red','DisplayName', 'Trajectory Plane');
hold off
alpha(0.4)

hold on
% Satellite Position and Path
% Plotting the first iteration
p = plot3(x(1),y(1),z(1),'b', 'HandleVisibility','off');
m = scatter3(x(1),y(1),z(1),'filled','b', 'DisplayName', 'Satellite');

%Plotting the Axis
plot3([0 25000], [0 0], [0 0], 'k', "Linewidth", 0.5, 'HandleVisibility','off');
plot3([0 0], [0 1.5*10^4], [0 0], 'k', "Linewidth", 0.5, 'HandleVisibility','off');
plot3([0 0], [0 0], [0 1.5*10^4], 'k', "Linewidth", 0.5, 'HandleVisibility','off');
%Plotting the Keplerian Elements
%a
xcenter = (max(x) + min(x))/2;
ycenter = (max(y) + min(y))/2;
zcenter = (max(z) + min(z))/2;
xend = xcenter + a(1);
a_vec = plot3([xcenter, xend], [ycenter, ycenter], [zcenter zcenter], 'g', "Linewidth", 2, 'DisplayName','a');

%e
r = [x(1), y(1), z(1)];
v = [xdot(1), ydot(1), zdot(1)];
H = cross(r, v);
mag_r = norm(r);
unit_r = r/mag_r;
e = (1/mu)*cross(v, H) - unit_r;
e_vec = plot3([0, e(1)], [0 e(2)], [0, e(3)], 'cyan', "Linewidth", 1, 'DisplayName','e');

%i
x_axis = [1 0 0]; y_axis = [0 1 0]; z_axis = [0 0 1];
h_vec = plot3([0, H(1)], [0 H(2)], [0, H(3)], 'k', "Linewidth", 0.5, 'HandleVisibility', 'off');
radius = norm(H)/2; 
start_deg = pi/2;
end_deg = (pi/2) + acos(dot((H/norm(H)), z_axis));
d_deg = 0.1;

i_x = []; i_y = []; i_z = [];

for i=start_deg:d_deg:end_deg
    i_x = [i_x 0];
    i_y = [i_y radius*cos(i)];
    i_z = [i_z radius*sin(i)];
end


i_ang = plot3(i_x, i_y, i_z, 'y', "Linewidth", 1, 'DisplayName', 'i');
%O

%w

%f



lgd = legend;
lgd.Location = 'northeast';

pause(inf);
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
    
    % Updating the title
    title(sprintf('Two-Body Trajectory - Keplerian\nTime: %0.2f secs', t(k)),...
    'Interpreter','Latex');
    % Delay
    pause(0.001)
    % Saving the figure
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',(t(k+1)-t(k))/2500);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',(t(k+1)-t(k))/2500);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, mag_e, e_vec, inc, O, w, f] = coordinate_conversion(x, y, z, xdot, ydot, zdot, mu)

a = []; mag_e = []; e_vec = []; inc = []; O = []; w = []; f = [];

unit_x = [1 0 0]; unit_y = [0 1 0]; unit_z = [0 0 1];

for k = 1:length(x)
    r = [x(k) y(k) z(k)];
    v = [xdot(k) ydot(k) zdot(k)];
    mag_r = norm(r);
    unit_r = r/mag_r;
    E = (0.5)*dot(v,v) - (mu/mag_r);
    a = [a, -(mu/(2*E))];
    H = cross(r, v);
    unit_h = H/(norm(H));
    e = (1/mu)*cross(v, H) - unit_r;
    norm_e = norm(e);
    mag_e = [mag_e norm_e];
    unit_e = e/(norm_e);
    inc = [inc acos(dot(unit_h, unit_z))];
    O = [O atan2(dot(unit_h, unit_x), dot(-unit_h, unit_y))];
    w = [w abs(atan2(dot(unit_e, unit_z), dot(cross(unit_h, unit_e), unit_z)))];

    if dot(r, v) < 0
        f = [f 2*pi - acos(dot(unit_r, unit_e))];
    else
        f = [f acos(dot(unit_r, unit_e))];
    end
       
end    
end 