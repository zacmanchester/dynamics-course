% Simple fixed-wing airplane simulation with LQR trajectory tracking.
% The airplane takes off, flies a left-hand racetrack and lands back on the
% runway.
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% where q = [scalar; vector] and maps body coordinates to world coordinates.
% Input:
%   u = [T; delta_a; delta_e; delta_r]
%
% Frames: inertial z up; body x forward, y toward the left wing, z up.
% Sign conventions: delta_e > 0 pitches the nose up, delta_a > 0 rolls the
% right wing down, delta_r > 0 yaws the nose to the left (toward +y).

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
p.Cm0 = -0.0135;  p.Cma  =  2.74;    p.Cmq  = -38.21;  p.Cmde = -0.99;
p.Cnb =  0.073;   p.Cnp  =  0.069;  p.Cnr  = -0.095;  p.Cnda =  0.011;  p.Cndr = 0.069;

% Actuator limits and runway
p.Tmax = 120;                % maximum thrust [N]
p.dmax = deg2rad(30);        % control-surface deflection limit
p.mu   = 0.4;                % runway friction (includes braking)

%% Flight plan
p.Vc       = 30;             % cruise speed [m/s]
p.hc       = 20;             % cruise altitude [m]
p.a_runway = 5;              % take-off acceleration [m/s^2]
p.a_brake  = 3.8;            % braking deceleration [m/s^2], just below mu*g
p.gamma_up   = deg2rad(8);   % climb angle
p.gamma_down = deg2rad(5);   % descent angle (shallow, so thrust stays positive)
p.bank     = deg2rad(45);    % bank angle in the turns
p.downwind = 310;            % downwind leg length [m]

Nx  = 13;           % quaternion state dimension
Nxr = 12;           % reduced/linearized state dimension
Nu  = 4;

%% Trim: straight and level flight at cruise speed
[x0,utrim,theta0] = trim_airplane(p);
p.theta0 = theta0;           % trim pitch angle, used by the reference trajectory
q0 = x0(4:7);

%% Linearize discrete-time RK4 dynamics about trim
fdisc = @(x,u) airplane_dynamics_rk4(x,u,p);
A = numerical_jacobian(@(x) fdisc(x,utrim), x0);
B = numerical_jacobian(@(u) fdisc(x0,u), utrim);

%% Reduced 12-state linearization
E0 = AirplaneMath.E(q0);
Ar = E0' * A * E0;
Br = E0' * B;

%% LQR controller
% Bryson's rule: weight = 1/(acceptable deviation)^2.
Qcost = diag(1./[5 5 5, 0.2 0.2 0.2, 3 3 3, 0.5 0.5 0.5].^2);
Rcost = diag(1./[40, 0.3 0.3 0.3].^2);
K = dlqr(Ar,Br,Qcost,Rcost);

umin = [0;      -p.dmax*ones(3,1)];
umax = [p.Tmax;  p.dmax*ones(3,1)];
controller = @(x,xref) feedback_controller(x,xref,utrim,K,umin,umax);

%% Simulation
[~,tphase] = reference_trajectory(0,p);
Tfinal = tphase(end) + 2.0;
Nt = round(Tfinal/p.h) + 1;
thist = (0:Nt-1)*p.h;

uhist    = zeros(Nu,Nt);
xhist    = zeros(Nx,Nt);
xrefhist = zeros(Nx,Nt);

% At rest on the runway, level, pointing along +x.
xhist(:,1) = [zeros(3,1); 1; 0; 0; 0; zeros(6,1)];

for k = 1:Nt-1
    xrefhist(:,k) = reference_trajectory(thist(k),p);
    uhist(:,k) = controller(xhist(:,k),xrefhist(:,k));
    xhist(:,k+1) = fdisc(xhist(:,k),uhist(:,k));

    % Stop if the airplane hits the ground away from the runway.
    if xhist(3,k+1) <= 0 && abs(xhist(2,k+1)) > 10
        fprintf('Crash at t = %.1f s\n',thist(k+1));
        Nt = k+1;
        thist    = thist(1:Nt);
        xhist    = xhist(:,1:Nt);
        uhist    = uhist(:,1:Nt);
        xrefhist = xrefhist(:,1:Nt);
        break
    end
end
xrefhist(:,Nt) = reference_trajectory(thist(Nt),p);
uhist(:,Nt) = controller(xhist(:,Nt),xrefhist(:,Nt));

%% Plots
figure('Name','Airplane position');
plot(thist,xhist(1:3,:),'LineWidth',1.5);
hold on;
plot(thist,xrefhist(1:3,:),'k--');
grid on;
xlabel('Time [s]');
ylabel('Position [m]');
legend('x','y','z','reference','Location','best');
title('Position');

figure('Name','Airplane controls');
subplot(2,1,1);
plot(thist,uhist(1,:),'LineWidth',1.2);
grid on;
ylabel('Thrust [N]');
title('Control inputs');
subplot(2,1,2);
plot(thist,rad2deg(uhist(2:4,:)),'LineWidth',1.2);
grid on;
xlabel('Time [s]');
ylabel('Deflection [deg]');
legend('\delta_a','\delta_e','\delta_r','Location','best');

figure('Name','Airplane air data');
Va    = vecnorm(xhist(8:10,:));
alpha = atan2(-xhist(10,:),xhist(8,:));
beta  = atan2(xhist(9,:),sqrt(xhist(8,:).^2 + xhist(10,:).^2));
subplot(2,1,1);
plot(thist,Va,'LineWidth',1.2);
grid on;
ylabel('Airspeed [m/s]');
title('Air data');
subplot(2,1,2);
plot(thist,rad2deg([alpha; beta]),'LineWidth',1.2);
grid on;
xlabel('Time [s]');
ylabel('Angle [deg]');
legend('\alpha','\beta','Location','best');

figure('Name','Airplane lift curve');
alphaGrid = deg2rad(-5:0.1:25);
CLlevel = 2*p.m*p.g/(p.rho*p.Vc^2*p.S);      % lift coefficient needed in level flight
CLturn  = CLlevel/cos(p.bank);               % ... and in the banked turn
plot(rad2deg(alphaGrid),lift_coefficient(alphaGrid,p),'LineWidth',1.5);
hold on;
plot(rad2deg(alpha),lift_coefficient(alpha,p),'r.');
yline(CLlevel,'k--','level flight');
yline(CLturn,'k--','turn');
grid on;
xlabel('\alpha [deg]');
ylabel('C_L');
legend('C_L(\alpha)','flown','Location','northwest');
title(sprintf('Lift curve at V = %.0f m/s',p.Vc));

figure('Name','Airplane trajectory');
plot3(xhist(1,:),xhist(2,:),xhist(3,:),'LineWidth',1.5);
hold on;
plot3(xrefhist(1,:),xrefhist(2,:),xrefhist(3,:),'k--');
grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
legend('flown','reference','Location','best');
title('3D trajectory');

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
