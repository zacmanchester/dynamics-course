% Fixed-wing airplane without any controller: the thrust and control
% surfaces are held at their trim values and the airplane is left on its own
% after a small pitch perturbation. Whether the motion dies out or grows is
% decided by the eigenvalues of the linearized dynamics A.
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% where q = [scalar; vector] and maps body coordinates to world coordinates.
% Input:
%   u = [T; delta_a; delta_e; delta_r]
%
% Frames: inertial z up; body x forward, y toward the left wing, z up.

clear; clc; close all;

%% Airplane parameters (small UAV)
p.m   = 11.0;
p.J   = diag([0.824, 1.135, 1.759]);
p.g   = 9.81;
p.rho = 1.268;
p.S   = 0.55;           % wing area
p.b   = 2.9;            % wingspan
p.c   = 0.19;           % mean aerodynamic chord
p.e   = 0.9;            % Oswald efficiency factor
p.AR  = p.b^2/p.S;      % aspect ratio
p.h   = 0.02;           % 50 Hz

% Aerodynamic coefficients (stability derivatives, z-up / y-left body frame)
p.CL0 =  0.23;    p.CLa  =  5.61;   p.CLq  = -7.95;   p.CLde = -0.13;
p.CD0 =  0.043;
p.alpha_s = deg2rad(10);     % stall angle
p.Ms      = 50;              % sharpness of the stall transition
p.CYb = -0.83;    p.CYda = -0.075;  p.CYdr = -0.19;
p.Clb =  0.13;    p.Clp  = -0.51;   p.Clr  = -0.25;   p.Clda =  0.17;   p.Cldr = 0.0024;
p.Cm0 = -0.0135;  p.Cma  =  -2.74;    p.Cmq  = -38.21;  p.Cmde = -0.99;
p.Cnb =  0.073;   p.Cnp  =  0.069;  p.Cnr  = -0.095;  p.Cnda =  0.011;  p.Cndr = 0.069;

% Actuator limits and runway
p.Tmax = 120;                % maximum thrust [N]
p.dmax = deg2rad(30);        % control-surface deflection limit
p.mu   = 0.4;                % runway friction (includes braking)

%% Flight condition
p.Vc = 30;                   % airspeed [m/s]
p.hc = 20;                   % altitude [m]

Nx  = 13;           % quaternion state dimension
Nxr = 12;           % reduced/linearized state dimension
Tfinal = 10.0;
Nt = round(Tfinal/p.h) + 1;
thist = (0:Nt-1)*p.h;

%% Trim: straight and level flight at cruise speed
[x0,utrim,theta0] = trim_airplane(p);
q0 = x0(4:7);

%% Linearize discrete-time RK4 dynamics about trim
fdisc = @(x,u) airplane_dynamics_rk4(x,u,p);
A = numerical_jacobian(@(x) fdisc(x,utrim), x0);

%% Reduced 12-state linearization
E0 = AirplaneMath.E(q0);
Ar = E0' * A * E0;

%% Eigenvalues of A
% The motion about trim is stable if all eigenvalues lie inside the unit circle.
figure('Name','Airplane eigenvalues');
phi = linspace(0,2*pi,200);
plot(cos(phi),sin(phi),'k--');
hold on; grid on; axis equal;
plot(real(eig(Ar)),imag(eig(Ar)),'bx','MarkerSize',8,'LineWidth',1.5);
xlabel('Re'); ylabel('Im');
legend('unit circle','eigenvalues of A','Location','best');
title(sprintf('Largest |eigenvalue| = %.4f',max(abs(eig(Ar)))));

%% Simulation: trim input held, 2 deg nose-up perturbation
xhist    = zeros(Nx,Nt);
xrefhist = zeros(Nx,Nt);          % undisturbed trim flight, for the animation

qpert = [cos(deg2rad(1)); 0; -sin(deg2rad(1)); 0];
xhist(:,1) = [x0(1:3); AirplaneMath.Lquat(q0)*qpert; x0(8:13)];

for k = 1:Nt-1
    xrefhist(:,k) = [p.Vc*thist(k); 0; p.hc; x0(4:13)];
    xhist(:,k+1) = fdisc(xhist(:,k),utrim);
end
xrefhist(:,Nt) = [p.Vc*thist(Nt); 0; p.hc; x0(4:13)];

%% Animation
animate_airplane(xhist,xrefhist,thist,p);


function [x0,u0,theta0] = trim_airplane(p)
    % Straight-and-level flight at speed p.Vc: find the pitch angle, elevator
    % and thrust for which the body accelerations vanish (Newton's method).
    z = zeros(3,1);                      % z = [theta; delta_e; T]
    for iter = 1:20
        res = trim_residual(z,p);
        Jz  = numerical_jacobian(@(z) trim_residual(z,p), z);
        z   = z - Jz\res;
    end
    theta0 = z(1);
    x0 = trim_state(theta0,p);
    u0 = [z(3); 0; z(2); 0];
end

function res = trim_residual(z,p)
    x  = trim_state(z(1),p);
    u  = [z(3); 0; z(2); 0];
    dx = airplane_dynamics(x,u,p);
    res = [dx(8); dx(10); dx(12)];       % forward, vertical and pitch accelerations
end

function x = trim_state(theta,p)
    % Level flight at pitch angle theta. The velocity is horizontal, so the
    % angle of attack equals theta.
    q = [cos(theta/2); 0; -sin(theta/2); 0];
    v = p.Vc*[cos(theta); 0; -sin(theta)];
    x = [0; 0; p.hc; q; v; zeros(3,1)];
end

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
