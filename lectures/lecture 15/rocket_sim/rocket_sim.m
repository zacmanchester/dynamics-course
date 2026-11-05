% Simple thrust-vectored rocket landing simulation with LQR tracking.
% The rocket starts tilted, off to the side of the pad and high up, and has to
% land upright on the pad with zero velocity.
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% where q = [scalar; vector] and maps body coordinates to world coordinates.
% Input:
%   u = [T; alpha; beta]  (thrust and the two gimbal angles)
%
% Frames: inertial z up; body z along the rocket axis (nose up).

clear; clc; close all;

%% Rocket parameters
p.m   = 1000;                     % mass [kg]
p.ell = 5;                        % center of mass to engine [m]
p.J   = diag([8300, 8300, 100]);  % slender body: small inertia about its axis
p.g   = 9.81;
p.h   = 0.02;                     % 50 Hz

p.Tmax = 2*p.m*p.g;               % maximum thrust [N]
p.gmax = deg2rad(10);             % gimbal angle limit

%% Landing plan
p.h0    = 60;                     % initial altitude [m]
p.Tland = 18;                     % time to reach the pad [s]

Nx  = 13;           % quaternion state dimension
Nxr = 12;           % reduced/linearized state dimension
Nu  = 3;
Tfinal = p.Tland + 3.0;
Nt = round(Tfinal/p.h) + 1;
thist = (0:Nt-1)*p.h;

%% Hover state and input
uhover = [p.m*p.g; 0; 0];
r0 = [0; 0; p.h0];
q0 = [1; 0; 0; 0];
v0 = zeros(3,1);
w0 = zeros(3,1);
x0 = [r0; q0; v0; w0];

%% Linearize discrete-time RK4 dynamics about hover
fdisc = @(x,u) rocket_dynamics_rk4(x,u,p);
A = numerical_jacobian(@(x) fdisc(x,uhover), x0);
B = numerical_jacobian(@(u) fdisc(x0,u), uhover);

%% Reduced 12-state linearization
E0 = RocketMath.E(q0);
Ar = E0' * A * E0;
Br = E0' * B;

%% LQR controller
% The gimbal cannot produce a torque about the rocket axis, so the roll angle
% and roll rate about that axis (reduced states 6 and 12) are left out.
idx = [1:5, 7:11];
% Bryson's rule: weight = 1/(acceptable deviation)^2.
Qcost = diag(1./[4 4 3, 0.1 0.1, 2 2 2, 0.2 0.2].^2);
Rcost = diag(1./[0.5*p.m*p.g, 0.03 0.03].^2);
K = zeros(Nu,Nxr);
K(:,idx) = dlqr(Ar(idx,idx),Br(idx,:),Qcost,Rcost);

umin = [0;      -p.gmax; -p.gmax];
umax = [p.Tmax;  p.gmax;  p.gmax];
controller = @(x,xref) feedback_controller(x,xref,uhover,K,umin,umax);

%% Simulation
uhist    = zeros(Nu,Nt);
xhist    = zeros(Nx,Nt);
xrefhist = zeros(Nx,Nt);

% Initial condition: off the pad in x and y, tilted 10 deg about a diagonal
% horizontal axis (so both gimbal angles are needed), descending at 3 m/s.
tilt = deg2rad(10);
tiltAxis = [1; 1; 0]/sqrt(2);
xhist(:,1) = [15; -10; p.h0; ...
              cos(tilt/2); sin(tilt/2)*tiltAxis; ...
              0; 0; -3; ...
              w0];

for k = 1:Nt-1
    xrefhist(:,k) = reference_trajectory(thist(k),p);
    uhist(:,k) = controller(xhist(:,k),xrefhist(:,k));
    if xhist(3,k) <= p.ell      % touched down: engine off
        uhist(:,k) = zeros(Nu,1);
    end
    xhist(:,k+1) = fdisc(xhist(:,k),uhist(:,k));
end
xrefhist(:,Nt) = reference_trajectory(thist(Nt),p);
uhist(:,Nt) = zeros(Nu,1);

%% Plots
figure('Name','Rocket position');
plot(thist,xhist(1:3,:),'LineWidth',1.5);
hold on;
plot(thist,xrefhist(1:3,:),'k--');
grid on;
xlabel('Time [s]');
ylabel('Position [m]');
legend('x','y','z','reference','Location','best');
title('Position');

figure('Name','Rocket controls');
subplot(2,1,1);
plot(thist,uhist(1,:)/(p.m*p.g),'LineWidth',1.2);
grid on;
ylabel('Thrust / weight');
title('Control inputs');
subplot(2,1,2);
plot(thist,rad2deg(uhist(2:3,:)),'LineWidth',1.2);
grid on;
xlabel('Time [s]');
ylabel('Gimbal angle [deg]');
legend('\alpha','\beta','Location','best');

figure('Name','Rocket trajectory');
plot3(xhist(1,:),xhist(2,:),xhist(3,:),'LineWidth',1.5);
hold on;
plot3(xrefhist(1,:),xrefhist(2,:),xrefhist(3,:),'k--');
grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
legend('flown','reference','Location','best');
title('3D trajectory');

%% Animation
animate_rocket(xhist,uhist,thist,p);


function A = numerical_jacobian(fun,x)
    % Central-difference Jacobian.
    fx = fun(x);
    m = numel(fx);
    n = numel(x);
    A = zeros(m,n);

    baseStep = eps^(1/3);
    for i = 1:n
        hi = baseStep*max(1,abs(x(i)));
        xp = x;
        xm = x;
        xp(i) = xp(i) + hi;
        xm(i) = xm(i) - hi;
        A(:,i) = (fun(xp)-fun(xm))/(2*hi);
    end
end
