% Flying car demo 2: the inverted ramp.
% The car accelerates on the ground and drives up a loop that turns from
% horizontal to vertical, climbs a vertical wall and curves back until it is
% upside down, held against the track by its suction fans. It crawls along
% the inverted straight at the top, drops off the end of the track, rolls
% 180 degrees about its longitudinal axis on its thrusters, lands on its
% wheels on the ground below and keeps driving.
%
% The script runs in three modes:
%   'drive'  : driving force and suction, thrusters off (no driving force
%              on the inverted straight: the car coasts along it)
%   'fly'    : LQR tracking of the flight reference (half roll in free fall,
%              brake the fall, lower)
%   'cruise' : all inputs off, the car rolls on
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% Input:
%   u = [T1; T2; T3; T4; T5; S; f]

clear; clc; close all;

%% Car parameters and loop
p = flyingcar_params();
p.loop.on        = true;
p.loop.x0        = 17;                % the loop starts here [m]
p.loop.R         = 6;                 % radius of the two quarter circles [m]
p.loop.Lv        = 4;                 % height of the vertical wall [m]
p.loop.Lt        = 6;                 % length of the inverted straight at the top [m]

%% Stunt plan
p.fdrive  = 0.10*p.m*p.g;         % driving force command per wheel [N]
p.S       = 1.3*p.m*p.g;          % suction while on the track: presses the car onto the ceiling [N]
p.T_roll  = 0.8;                  % duration of the half roll [s]
p.a_brake = 0.35*p.g;             % upward acceleration to stop the fall
p.vdes    = 1.5;                  % descent speed when lowering [m/s]

Nx  = 13;           % quaternion state dimension
Nxr = 12;           % reduced/linearized state dimension
Nu  = 7;
Tfinal = 17.0;
Nt = round(Tfinal/p.h) + 1;
thist = (0:Nt-1)*p.h;

%% Hover state and input (in the air, no tire forces)
uhover = [p.m*p.g/4*ones(4,1); 0; 0; 0];
r0 = [0; 0; 3];
q0 = [1; 0; 0; 0];
v0 = zeros(3,1);
w0 = zeros(3,1);
x0 = [r0; q0; v0; w0];

%% Linearize discrete-time RK4 dynamics about hover
fdisc = @(x,u) flyingcar_dynamics_rk4(x,u,p);
A = numerical_jacobian(@(x) fdisc(x,uhover), x0);
B = numerical_jacobian(@(u) fdisc(x0,u), uhover);

%% Reduced 12-state linearization
E0 = FlyingCarMath.E(q0);
Ar = E0' * A * E0;
Br = E0' * B;

%% LQR controller
% Only the five thrusters are used in the air; suction and driving force
% stay zero (their rows of K are left at zero). The thrusters cannot produce
% a yaw torque, so the yaw angle and yaw rate (reduced states 6 and 12) are
% left out.
iu = 1:5;
ix = [1:5, 7:11];
% Bryson's rule: weight = 1/(acceptable deviation)^2.
Qcost = diag(1./[0.5 0.5 0.5, 0.05 0.05, 1 1 1, 0.2 0.2].^2);
Rcost = diag(1./(0.3*p.m*p.g*ones(1,5)).^2);
K = zeros(Nu,Nxr);
K(iu,ix) = dlqr(Ar(ix,ix),Br(ix,iu),Qcost,Rcost);

umin = zeros(Nu,1);
umax = [p.Tmax*ones(4,1); p.T5max; p.Smax; inf];
controller = @(x,xref,uref,K) feedback_controller(x,xref,uref,K,umin,umax);

% While rolling, the position feedback is kept only along the thrust axis
% (reduced states 3 and 9): a position error sideways would otherwise be
% fought by rolling the car, against the roll we are asking for. It is also
% softened, so that height corrections do not saturate the thrusters and
% spoil the roll torque.
Kroll = K;
Kroll(:,[1 2 7 8]) = 0;
Kroll(:,[3 9]) = 0.25*Kroll(:,[3 9]);

%% Simulation
uhist    = zeros(Nu,Nt);
xhist    = zeros(Nx,Nt);
xrefhist = nan(Nx,Nt);
Nhist    = zeros(4,Nt);           % tire normal forces
mode     = 'drive';

% Initial condition: at rest on the ground, before the loop.
xhist(:,1) = [0; 0; p.hw; q0; v0; w0];

for k = 1:Nt-1
    x = xhist(:,k);
    [~,~,Nhist(:,k)] = tire_forces(x,0,p);
    onGround = any(Nhist(:,k) > 0);

    switch mode
        case 'drive'
            onCeiling = x(3) > p.loop.R + p.loop.Lv && x(1) < p.loop.x0;
            pastEnd   = onCeiling && x(1) - p.a < p.loop.x0 - p.loop.Lt;   % front wheels off the track
            % Throttle off on the ceiling (the car coasts along it), suction
            % off once the front wheels leave the track (the fans lose their seal).
            uhist(:,k) = [0; 0; 0; 0; 0; p.S*~pastEnd; p.fdrive*~onCeiling];
            if ~onGround && onCeiling                     % dropped off the end of the track
                mode = 'fly';
                t_launch = thist(k);
                r_launch = x(1:3);
                v_launch = FlyingCarMath.qtoQ(x(4:7))*x(8:10);   % world frame
            end
        case 'fly'
            [xrefhist(:,k),uref,rolling] = fly_reference(thist(k) - t_launch,r_launch,v_launch,p);
            if rolling
                uhist(:,k) = controller(x,xrefhist(:,k),uref,Kroll);
            else
                uhist(:,k) = controller(x,xrefhist(:,k),uref,K);
            end
            if onGround                           % landed
                mode = 'cruise';
            end
        case 'cruise'
            uhist(:,k) = zeros(Nu,1);
    end

    xhist(:,k+1) = fdisc(x,uhist(:,k));
end
[~,~,Nhist(:,Nt)] = tire_forces(xhist(:,Nt),0,p);

%% Plots
figure('Name','Ramp position');
plot(thist,xhist(1:3,:),'LineWidth',1.5);
hold on;
plot(thist,xrefhist(1:3,:),'k--');
grid on;
xlabel('Time [s]');
ylabel('Position [m]');
legend('x','y','z','reference','Location','best');
title('Position');

figure('Name','Ramp attitude');
Qz = zeros(3,Nt);                  % body z axis in the world frame
for k = 1:Nt
    Q = FlyingCarMath.qtoQ(xhist(4:7,k));
    Qz(:,k) = Q(:,3);
end
plot(thist,Qz,'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('Body z axis in world frame');
legend('x component','y component','z component (1 = upright, -1 = inverted)','Location','best');
title('Attitude');

figure('Name','Ramp controls');
subplot(2,1,1);
plot(thist,uhist(1:5,:)/(p.m*p.g),'LineWidth',1.2);
grid on;
ylabel('Thrust / weight');
legend('T_1','T_2','T_3','T_4','T_5','Location','best');
title('Thruster commands');
subplot(2,1,2);
plot(thist,Nhist/(p.m*p.g),'LineWidth',1.2);
grid on;
xlabel('Time [s]');
ylabel('Tire normal force / weight');
legend('N_1','N_2','N_3','N_4','Location','best');

figure('Name','Ramp side view');
[xl,zl] = loop_profile(p);
plot(xl,zl,'k-','LineWidth',2);
hold on;
plot(xhist(1,:),xhist(3,:),'LineWidth',1.5);
grid on; axis equal;
xlabel('x [m]'); ylabel('z [m]');
legend('track','center of mass','Location','best');
title('Side view');

%% Animation
animate_flyingcar(xhist,uhist,thist,p);


function [xref,uref,rolling] = fly_reference(tau,r_launch,v_launch,p)
    % Reference state and feedforward input tau seconds after leaving the
    % track. Height: free fall while rolling, then brake the fall at a_brake
    % until the descent speed is vdes, then lower at vdes to the ground.
    % Attitude: the launch pitch (nose along the track end), plus a half
    % roll about body x that turns the car upright.
    z1  = r_launch(3) + v_launch(3)*p.T_roll - 0.5*p.g*p.T_roll^2;   % end of free fall
    vz1 = v_launch(3) - p.g*p.T_roll;
    t_brake = (-p.vdes - vz1)/p.a_brake;
    z2  = z1 + vz1*t_brake + 0.5*p.a_brake*t_brake^2;                % end of braking

    rolling = tau < p.T_roll;
    if rolling                                          % free fall, rolling
        z    = r_launch(3) + v_launch(3)*tau - 0.5*p.g*tau^2;
        zdot = v_launch(3) - p.g*tau;
        T0   = 0;
    elseif tau < p.T_roll + t_brake                     % braking the fall
        s    = tau - p.T_roll;
        z    = z1 + vz1*s + 0.5*p.a_brake*s^2;
        zdot = vz1 + p.a_brake*s;
        T0   = p.m*(p.g + p.a_brake)/4;
    else                                                % lowering
        s    = tau - p.T_roll - t_brake;
        z    = max(z2 - p.vdes*s, p.hw - 0.05);     % a bit below the resting
        zdot = -p.vdes*(z > p.hw - 0.05);          % height: wheels take the weight
        T0   = p.m*p.g/4;
    end

    % Half roll about body x, from inverted to upright.
    s       = tau/p.T_roll;
    phi     = pi*smoothstep(s);
    phidot  = pi*dsmoothstep(s)/p.T_roll;
    phiddot = pi*ddsmoothstep(s)/p.T_roll^2;

    th = pi;                                            % upside down, nose along -x
    qpitch = [cos(th/2); 0; -sin(th/2); 0];
    qroll  = [cos(phi/2); sin(phi/2); 0; 0];
    qref   = FlyingCarMath.Lquat(qpitch)*qroll;
    Qref   = FlyingCarMath.qtoQ(qref);

    xref = [r_launch(1) + v_launch(1)*tau; 0; z; ...
            qref; ...
            Qref'*[v_launch(1); 0; zdot]; ...
            phidot; 0; 0];

    % Feedforward: the thrusters carry the car only while they point up,
    % and the pair on one side fires harder to produce the roll torque
    % J_xx*phiddot = b*(T1 + T4 - T2 - T3).
    T0 = T0*max(Qref(3,3),0);
    dT = p.J(1,1)*phiddot/(2*p.b);
    uref = [T0 + max(dT,0); T0 + max(-dT,0); T0 + max(-dT,0); T0 + max(dT,0); 0; 0; 0];
end

function s = smoothstep(s)
    % 0 for s<=0, 1 for s>=1, smooth in between.
    s = min(max(s,0),1);
    s = 3*s^2 - 2*s^3;
end

function ds = dsmoothstep(s)
    % First derivative of smoothstep.
    s = min(max(s,0),1);
    ds = 6*s - 6*s^2;
end

function dds = ddsmoothstep(s)
    % Second derivative of smoothstep.
    dds = (6 - 12*s)*(s > 0 && s < 1);
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
