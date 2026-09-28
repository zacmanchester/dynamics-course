clear all;
close all;

%We're going to simulate our ant on a frictionless record (you can think of
%it as made of ice).

c = 0.0; %friction coefficient

%Assume record B frame and inertial N frame are aligned at t=0
%Assume record is spinning at omega = 33.333 RPM = 3.5 rad/s about z

omega = [0; 0; 3.5];

%Assume ant starts out stationary w.r.t. the record at r = [0.01; 0; 0] m

r0 = [.01; 0; 0];

rdot0_B = [0; 0; 0];

rdot0_N = hat(omega)*r0;

x0_N = [r0; rdot0_N];
x0_B = [r0; rdot0_B];

Tfinal = 10;

[t_N, x_N] = ode45(@(t,x)ant_inertial(t,x,omega,c),[0 Tfinal], x0_N);

[t_B, x_B] = ode45(@(t,x)ant_record(t,x,omega,c), [0 Tfinal], x0_B);


figure(1)
plot(x_N(:,1),x_N(:,2), 'LineWidth', 2);
hold on;
viscircles([0, 0], 0.3, 'color', 'k');
axis equal;


figure(2);
plot(x_B(:,1),x_B(:,2), 'LineWidth', 2);
hold on;
viscircles([0, 0], 0.3, 'color', 'k');
axis equal;