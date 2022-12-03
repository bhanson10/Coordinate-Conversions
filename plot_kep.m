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

mu = 3.986*10^5; % Earth's gravitational parameter [km^3/s^2]
[a, mag_e, inc, O, w, f] = coordinate_conversion(x, y, z, xdot, ydot, zdot, mu);

tiledlayout(3,1) 

% Top plot
nexttile
hold on; grid on;
a_plot = plot(t(1),a(1), 'm', 'LineWidth', 2, 'HandleVisibility','off');
a_cur_plot = scatter(t(1),a(1), 'm', 'filled', 'DisplayName', 'Semimajor Axis');
title('Keplerian Orbital Elements over Time', 'Interpreter','Latex', 'FontSize', 18);
ylim([1.402*10^4, 1.422*10^4]);
xlim([0, t(length(t))]);
lgd = legend;
lgd.Location = "southeast";
lgd.FontSize = 10;
xlabel('t (s)', 'Interpreter', 'Latex', 'FontSize', 16)
ylabel('Length (km)', 'Interpreter', 'Latex', 'FontSize', 16)

% Middle plot
nexttile
hold on; grid on;
e_plot = plot(t(1),mag_e(1), 'c', 'LineWidth', 2, 'HandleVisibility','off');
e_cur_plot = scatter(t(1),mag_e(1), 'filled', 'c', 'DisplayName', 'Eccentricity');
ylim([.411, .42]);
xlim([0, t(length(t))]);
lgd = legend;
lgd.Location = "southeast";
lgd.FontSize = 10;
xlabel('t (s)', 'Interpreter', 'Latex', 'FontSize', 16)

% Bottom plot
nexttile
hold on; grid on;
inc_plot = plot(t(1),inc(1), 'b', 'Linewidth', 2, 'HandleVisibility', 'off');
inc_cur_plot = scatter(t(1),inc(1), 'filled','b', 'DisplayName', 'Inclination'); 
O_plot = plot(t(1),O(1), 'r', 'Linewidth', 2, 'HandleVisibility', 'off'); 
O_cur_plot = scatter(t(1),O(1), 'r', 'filled', 'DisplayName', 'Longitude of Asc. Node');
w_plot = plot(t(1),w(1), 'g', 'Linewidth', 2, 'HandleVisibility','off');
w_cur_plot = scatter(t(1),w(1), 'g', 'filled', 'DisplayName', 'Arg. of Periapsis');
f_plot = plot(t(1),f(1), 'k', 'Linewidth', 2, 'HandleVisibility', 'off');
f_cur_plot = scatter(t(1),f(1), 'k', 'filled', 'DisplayName', 'True Anomaly');
ylim([-2*pi*1.1, 2*pi*1.1]);
xlim([0, t(length(t))]);
lgd = legend;
lgd.Location = "southeast";
lgd.FontSize = 6;
xlabel('t (s)', 'Interpreter', 'Latex', 'FontSize', 16)
ylabel('Angle (rad)', 'Interpreter', 'Latex', 'FontSize', 16)

% Create file name variable
filename = 'kep-elem.gif';

pause(5)
% Iterating through the length of the time array
for k = 1:length(t)-1
    
    % Updating the point
    a_cur_plot.XData = t(k); a_cur_plot.YData = a(k);
    e_cur_plot.XData = t(k); e_cur_plot.YData = mag_e(k);
    inc_cur_plot.XData = t(k); inc_cur_plot.YData = inc(k);
    O_cur_plot.XData = t(k); O_cur_plot.YData = O(k);
    w_cur_plot.XData = t(k); w_cur_plot.YData = w(k);
    f_cur_plot.XData = t(k); f_cur_plot.YData = f(k);

    % Updating the lines
    a_plot.XData = t(1:k); a_plot.YData = a(1:k);
    e_plot.XData = t(1:k); e_plot.YData = mag_e(1:k);
    inc_plot.XData = t(1:k); inc_plot.YData = inc(1:k);
    O_plot.XData = t(1:k); O_plot.YData = O(1:k);
    w_plot.XData = t(1:k); w_plot.YData = w(1:k);
    f_plot.XData = t(1:k); f_plot.YData = f(1:k);
    
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, mag_e, inc, O, w, f] = coordinate_conversion(x, y, z, xdot, ydot, zdot, mu)

a = []; mag_e = []; inc = []; O = []; w = []; f = [];

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