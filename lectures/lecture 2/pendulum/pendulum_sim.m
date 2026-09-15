clear all;
close all;

Tfinal = 60;

%initial conditions
x0_polar = [0; 0];
x0_cartesian = [1; 0; 0; 0];

[t_polar, x_polar] = ode45(@pendulum_polar,[0 Tfinal],x0_polar);

[t_cartesian, x_cartesian] = ode45(@pendulum_cartesian,[0 Tfinal],x0_cartesian);
%[t_cartesian, x_cartesian] = ode45(@pendulum_cartesian_pd,[0 Tfinal],x0_cartesian);

figure(1)
plot(t_polar, x_polar(:,1), LineWidth=2)

theta_cartesian = zeros(length(t_cartesian),1);
r_cartesian = zeros(length(t_cartesian),1);
for k = 1:length(t_cartesian)
    theta_cartesian(k) = atan2(x_cartesian(k,2),x_cartesian(k,1));
    r_cartesian(k) = norm(x_cartesian(k,1:2));
end

hold on
plot(t_cartesian, theta_cartesian, 'r', LineWidth=2)

figure(3)
plot(x_cartesian(:,1),x_cartesian(:,2))

figure(2)
plot(t_cartesian, r_cartesian)

