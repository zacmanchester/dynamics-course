clear all;
close all;

m = 1;
l = 1;
g = 9.8;

Tfinal = 60;

%initial conditions
x0 = [pi/4; 0];

[t_ode45, x_ode45] = ode45(@pendulum_dynamics,[0 Tfinal],x0);

figure(1)
plot(t_ode45, x_ode45(:,1), LineWidth=2)

dt = 0.001;
t_euler = 0:dt:Tfinal;

x_euler = zeros(2,length(t_euler));
x_euler(:,1) = x0;
for k = 1:(length(t_euler)-1)
    x_euler(:,k+1) = x_euler(:,k) + dt*pendulum_dynamics(t_euler(k),x_euler(:,k));
end

hold on
plot(t_euler, x_euler(1,:), 'r', LineWidth=2)

figure(2)
plot(x_ode45(:,1),x_ode45(:,2))

E0 = m*g*l*(1-cos(x0(1))) + 0.5*m*l^2*x0(2)^2;
theta_max = acos(1-E0/(m*g*l));
thetadot_max = sqrt(2*E0/(m*l^2));

hold on
rectangle('Position', [-theta_max, -thetadot_max, 2*theta_max, 2*thetadot_max], 'EdgeColor', 'g', 'LineWidth', 2);