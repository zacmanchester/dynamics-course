% Flying car demo 1: the kickflip.
% The car accelerates on the ground, jumps on its thrusters over a row of
% five parked cars, rolls a full 360 degrees about its longitudinal axis
% while flying over them, lands on its wheels and brakes to a stop.
%
% The script runs in three modes:
%   'drive' : driving force, suction and the rear thruster
%   'fly'   : LQR tracking of the jump reference: a gentle boost, the roll
%             in free flight (upside down over the tallest car), a gentle
%             braking of the fall
%   'brake' : braking force, thrusters off
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% Input:
%   u = [T1; T2; T3; T4; T5; S; f]

clear; clc; close all;

%% Car parameters
p = flyingcar_params();

%% Stunt plan
p.fdrive = 0.25*p.m*p.g;          % driving force command per wheel [N]
p.fbrake = 0.25*p.m*p.g;          % braking force per wheel [N]
p.T5     = 0.5*p.m*p.g;           % rear thruster during the run-up [N]
p.S      = 0.5*p.m*p.g;           % suction while on the ground [N]
p.x_jump = 3;                     % the jump starts here [m]
p.T_roll = 1.8;                   % duration of the roll, done in free flight [s]

% Five parked cars to jump over; the one in the middle is twice as tall.
p.cars.x = p.x_jump + 14 + (0:4)*4.5;     % centers [m]
p.cars.H = [1.3, 1.3, 2.6, 1.3, 1.3];     % heights [m]
p.cars.L = 4; p.cars.W = 1.8;

Nx  = 13;           % quaternion state dimension
Nxr = 12;           % reduced/linearized state dimension
Nu  = 7;
Tfinal = 10.0;
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

% Initial condition: at rest on the ground.
xhist(:,1) = [0; 0; p.hw; q0; v0; w0];

for k = 1:Nt-1
    x = xhist(:,k);
    [~,~,Nhist(:,k)] = tire_forces(x,0,p);
    onGround = any(Nhist(:,k) > 0);

    switch mode
        case 'drive'
            uhist(:,k) = [0; 0; 0; 0; p.T5; p.S; p.fdrive];
            if x(1) >= p.x_jump                   % jump!
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
            if onGround && thist(k) - t_launch > 1  % landed
                mode = 'brake';
            end
        case 'brake'
            uhist(:,k) = [0; 0; 0; 0; 0; p.S; -p.fbrake*(x(8) > 0.1)];
    end

    xhist(:,k+1) = fdisc(x,uhist(:,k));
end
[~,~,Nhist(:,Nt)] = tire_forces(xhist(:,Nt),0,p);

%% Plots
figure('Name','Kickflip position');
plot(thist,xhist(1:3,:),'LineWidth',1.5);
hold on;
plot(thist,xrefhist(1:3,:),'k--');
grid on;
xlabel('Time [s]');
ylabel('Position [m]');
legend('x','y','z','reference','Location','best');
title('Position');

figure('Name','Kickflip roll angle');
roll    = 2*atan2(xhist(5,:),xhist(4,:));
rollref = 2*atan2(xrefhist(5,:),xrefhist(4,:));
plot(thist,rad2deg(unwrap(roll)),'LineWidth',1.5);
hold on;
plot(thist,rad2deg(unwrap(rollref)),'k--');
grid on;
xlabel('Time [s]');
ylabel('Roll angle [deg]');
legend('flown','reference','Location','best');
title('Roll angle');

figure('Name','Kickflip controls');
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

figure('Name','Kickflip side view');
plot(xhist(1,:),xhist(3,:),'LineWidth',1.5);
hold on;
for i = 1:numel(p.cars.x)
    rectangle('Position',[p.cars.x(i)-p.cars.L/2, 0, p.cars.L, p.cars.H(i)],'FaceColor',[0.6 0.6 0.6]);
end
grid on; axis equal;
xlabel('x [m]'); ylabel('z [m]');
title('Side view (center of mass)');

%% Animation
animate_flyingcar(xhist,uhist,thist,p);


function [xref,uref,rolling] = fly_reference(tau,r_launch,v_launch,p)
    % Reference state and feedforward input tau seconds after the jump.
    % The car is thrown up gently on its thrusters, coasts through one full
    % roll in free flight (upside down right over the tallest car, at the
    % top of the arc), and is then caught gently by its thrusters again.
    % Forward speed is kept throughout.
    t_mid = (p.cars.x(3) - r_launch(1))/v_launch(1);    % passing over the tallest car
    t1 = t_mid - p.T_roll/2;        % boost ends, roll starts
    t2 = t_mid + p.T_roll/2;        % roll ends, braking starts
    t3 = t2 + t1;                   % braking ends, back at the launch height
    vz0 = p.g*p.T_roll/2;           % vertical speed to coast up and down during the roll
    a   = vz0/t1;                   % boost and braking acceleration
    z1  = r_launch(3) + 0.5*a*t1^2; % height at the start of the roll
    z_rest = p.hw - 0.05;           % a bit below the resting height: wheels take the weight

    rolling = false;
    if tau < t1                                         % boost
        z    = r_launch(3) + 0.5*a*tau^2;
        zdot = a*tau;
        T0   = p.m*(p.g + a)/4;
    elseif tau < t2                                     % free flight, rolling
        s    = tau - t1;
        z    = z1 + vz0*s - 0.5*p.g*s^2;
        zdot = vz0 - p.g*s;
        T0   = 0;
        rolling = true;
    elseif tau < t3                                     % braking the fall
        s    = tau - t2;
        z    = z1 - vz0*s + 0.5*a*s^2;
        zdot = -vz0 + a*s;
        T0   = p.m*(p.g + a)/4;
    else                                                % settle on the ground
        z    = z_rest;
        zdot = 0;
        T0   = p.m*p.g/4;
    end

    % Roll angle: one full turn about body x during the free flight.
    s       = (tau - t1)/p.T_roll;
    phi     = 2*pi*smoothstep(s);
    phidot  = 2*pi*dsmoothstep(s)/p.T_roll;
    phiddot = 2*pi*ddsmoothstep(s)/p.T_roll^2;

    qref = [cos(phi/2); sin(phi/2); 0; 0];
    Qref = FlyingCarMath.qtoQ(qref);

    xref = [r_launch(1) + v_launch(1)*tau; 0; z; ...
            qref; ...
            Qref'*[v_launch(1); 0; zdot]; ...
            phidot; 0; 0];

    % Feedforward: the thrusters carry the car while they point up, and the
    % pair on one side fires harder to produce the roll torque
    % J_xx*phiddot = b*(T1 + T4 - T2 - T3).
    T0 = T0*max(cos(phi),0);
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
