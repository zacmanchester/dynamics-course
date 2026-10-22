% Simple quadrotor simulation with LQR stabilization.
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% where q = [scalar; vector] and maps body coordinates to world coordinates.

clear; clc; close all;

%% Quadrotor parameters
p.m   = 0.5;
p.ell = 0.1750;
p.J   = diag([0.0023, 0.0023, 0.004]);
p.g   = 9.81;
p.kt  = 1.0;
p.km  = 0.0245;
p.h   = 0.05;       % 20 Hz

Nx  = 13;           % quaternion state dimension
Nxr = 12;           % reduced/linearized state dimension
Nu  = 4;
Tfinal = 5.0;
Nt = round(Tfinal/p.h) + 1;
thist = (0:Nt-1)*p.h;

%% Initial conditions and hover input
uhover = (p.m*p.g/4)*ones(4,1);
r0 = [0; 0; 1.0];
q0 = [1; 0; 0; 0];
v0 = zeros(3,1);
w0 = zeros(3,1);
x0 = [r0; q0; v0; w0];

%% Linearize discrete-time RK4 dynamics about hover
fdisc = @(x,u) quad_dynamics_rk4(x,u,p);
A = numerical_jacobian(@(x) fdisc(x,uhover), x0);
B = numerical_jacobian(@(u) fdisc(x0,u), uhover);

%% Reduced 12-state linearization
E0 = QuadrotorMath.E(q0);
Ar = E0' * A * E0;
Br = E0' * B;

%% LQR controller
Qcost = eye(Nxr);
Rcost = 0.1*eye(Nu);
K = dlqr(Ar,Br,Qcost,Rcost);

controller = @(x) feedback_controller(x,x0,r0,v0,w0,uhover,K);

%% Simulation
uhist = zeros(Nu,Nt);
xhist = zeros(Nx,Nt);

% perturbed initial condition:
%   r = r0 + randn(3)
%   attitude Rodrigues parameter phi = [1;0;0] + 0.5*randn(3)
xhist(:,1) = [r0 + randn(3,1); ...
              QuadrotorMath.Lquat(q0)*QuadrotorMath.rptoq([1;0;0] + 0.1*randn(3,1)); ...
              v0; w0];

for k = 1:Nt-1
    uhist(:,k) = controller(xhist(:,k));
    xhist(:,k+1) = fdisc(xhist(:,k),uhist(:,k));
end
uhist(:,Nt) = controller(xhist(:,Nt));

%% Plots
figure('Name','Quadrotor position');
plot(thist,xhist(1:3,:),'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('Position [m]');
legend('x','y','z','Location','best');
title('Position');

figure('Name','Quadrotor controls');
plot(thist,uhist,'LineWidth',1.2);
grid on;
xlabel('Time [s]');
ylabel('Rotor thrust command');
legend('u_1','u_2','u_3','u_4','Location','best');
title('Control inputs');

figure('Name','Quadrotor trajectory');
plot3(xhist(1,:),xhist(2,:),xhist(3,:),'LineWidth',1.5);
grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('3D trajectory');

%% Animation
animate_quadrotor(xhist,thist,p.ell);


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
