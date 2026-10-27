% QUADROTOR_KEYBOARD  Real-time keyboard-controlled quadrotor simulation.
%
% Controls (hold keys):
%   Up / Down arrows    : +x / -x world-frame reference velocity
%   Left / Right arrows : +y / -y world-frame reference velocity
%   a / s               : +z / -z world-frame reference velocity
%   d / f               : +yaw / -yaw reference angular velocity
%   escape or q          : quit
%
% Releasing a translation key commands zero velocity. The position reference
% is obtained by integrating the commanded velocity, so the quadrotor comes
% to a stop and holds its new location. Yaw works the same way: yaw rate is
% integrated into a yaw-angle reference, which is held when d/f are released.
%
% Requires the other files in this folder and Control System Toolbox (dlqr).

clear; clc; close all;

%% Quadrotor parameters
p.m   = 0.5;
p.ell = 0.1750;
p.J   = diag([0.0023, 0.0023, 0.004]);
p.g   = 9.81;
p.kt  = 1.0;
p.km  = 0.0245;
p.h   = 0.05;       % 20 Hz simulation/control update

Nx  = 13;
Nxr = 12;
Nu  = 4;

%% Keyboard command magnitudes
vxy_cmd  = 1.0;       % [m/s]
vz_cmd   = 0.6;       % [m/s]
yaw_cmd  = pi/3;      % [rad/s]

%% Hover state and input
r0 = [0; 0; 1.0];
q0 = [1; 0; 0; 0];
v0 = zeros(3,1);
w0 = zeros(3,1);
x0 = [r0; q0; v0; w0];
uhover = (p.m*p.g/4)*ones(4,1);

%% Linearize discrete-time RK4 dynamics about hover and form LQR gain
fdisc = @(x,u) quad_dynamics_rk4(x,u,p);
A = numerical_jacobian(@(x) fdisc(x,uhover), x0);
B = numerical_jacobian(@(u) fdisc(x0,u), uhover);

E0 = QuadrotorMath.E(q0);
Ar = E0' * A * E0;
Br = E0' * B;

Qcost = eye(Nxr);
Rcost = 0.1*eye(Nu);
K = dlqr(Ar,Br,Qcost,Rcost);

%% Initial state and moving reference
x = x0;
rref = r0;
yawref = 0;

%% Interactive figure
fig = figure('Name','Keyboard quadrotor control', ...
             'NumberTitle','off', ...
             'Color','w');
ax = axes(fig);
hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); axis(ax,'vis3d');
xlabel(ax,'x [m]'); ylabel(ax,'y [m]'); zlabel(ax,'z [m]');
view(ax,35,25);

set(fig,'WindowKeyPressFcn',@key_press);
set(fig,'WindowKeyReleaseFcn',@key_release);
setappdata(fig,'keys',struct());
setappdata(fig,'quit',false);

% Ground plane marker and trajectory.
plot3(ax,[-20 20],[0 0],[0 0],':','Color',[0.75 0.75 0.75]);
plot3(ax,[0 0],[-20 20],[0 0],':','Color',[0.75 0.75 0.75]);
traj = animatedline(ax,'LineWidth',1.2,'LineStyle','--');
refMarker = plot3(ax,rref(1),rref(2),rref(3),'rx','MarkerSize',10,'LineWidth',1.5);

% Body arms.
arm1 = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'k-','LineWidth',2.5);
arm2 = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'k-','LineWidth',2.5);

% Propeller circles. Four small blue circles in the body x-y plane.
Ncircle = 30;
theta = linspace(0,2*pi,Ncircle);
propRadius = 0.22*p.ell;
propHandles = gobjects(4,1);
for i = 1:4
    propHandles(i) = plot3(ax,nan(1,Ncircle),nan(1,Ncircle),nan(1,Ncircle), ...
                           'b-','LineWidth',1.5);
end

instructionText = sprintf([ ...
    'Hold keys:  \\uparrow/\\downarrow = +/-x   \\leftarrow/\\rightarrow = +/-y   ' ...
    'a/s = +/-z   d/f = +/-yaw   q or Esc = quit']);
title(ax,instructionText,'FontSize',10);

statusText = text(ax,0,0,0,'','FontName','FixedWidth','VerticalAlignment','bottom');

followRadiusXY = 2.5;
followBelowZ = 1.0;
followAboveZ = 2.0;

%% Real-time simulation loop
while ishandle(fig) && ~getappdata(fig,'quit')
    ticStep = tic;

    keys = getappdata(fig,'keys');

    % World-frame translational velocity reference.
    vrefN = zeros(3,1);
    if key_down(keys,'uparrow'),    vrefN(1) = vrefN(1) + vxy_cmd; end
    if key_down(keys,'downarrow'),  vrefN(1) = vrefN(1) - vxy_cmd; end
    if key_down(keys,'leftarrow'),  vrefN(2) = vrefN(2) + vxy_cmd; end
    if key_down(keys,'rightarrow'), vrefN(2) = vrefN(2) - vxy_cmd; end
    if key_down(keys,'a'),          vrefN(3) = vrefN(3) + vz_cmd;  end
    if key_down(keys,'s'),          vrefN(3) = vrefN(3) - vz_cmd;  end

    yawRateRef = 0;
    if key_down(keys,'d'), yawRateRef = yawRateRef + yaw_cmd; end
    if key_down(keys,'f'), yawRateRef = yawRateRef - yaw_cmd; end

    % Integrate velocity commands into a pose reference.
    rref = rref + p.h*vrefN;
    yawref = wrap_angle(yawref + p.h*yawRateRef);

    % Reference attitude is level with the commanded yaw.
    qref = [cos(0.5*yawref); 0; 0; sin(0.5*yawref)];
    Qref = QuadrotorMath.qtoQ(qref);

    % Express translational errors in the reference-yaw frame. This uses the
    % yaw symmetry of the hover linearization, while retaining the same K.
    dr = Qref'*(x(1:3) - rref);
    q = x(4:7)/norm(x(4:7));
    phi = QuadrotorMath.qtorp(QuadrotorMath.Lquat(qref)'*q);
    vrefB = Qref'*vrefN;
    dv = x(8:10) - vrefB;
    dw = x(11:13) - [0;0;yawRateRef];

    dxr = [dr; phi; dv; dw];
    u = uhover - K*dxr;

    % A rotor cannot produce negative thrust.
    u = max(u,0);

    % Advance one fixed simulation step.
    x = fdisc(x,u);

    %% Update graphics
    r = x(1:3);
    q = x(4:7)/norm(x(4:7));
    Q = QuadrotorMath.qtoQ(q);

    bodyPts1 = p.ell*[1, -1; 0, 0; 0, 0];
    bodyPts2 = p.ell*[0, 0; 1, -1; 0, 0];
    p1 = r + Q*bodyPts1;
    p2 = r + Q*bodyPts2;
    set(arm1,'XData',p1(1,:),'YData',p1(2,:),'ZData',p1(3,:));
    set(arm2,'XData',p2(1,:),'YData',p2(2,:),'ZData',p2(3,:));

    rotorCenters = p.ell*[ 1, -1,  0,  0; ...
                           0,  0,  1, -1; ...
                           0,  0,  0,  0];
    circleOffsets = propRadius*[cos(theta); sin(theta); zeros(1,Ncircle)];
    for i = 1:4
        propBody = rotorCenters(:,i) + circleOffsets;
        propWorld = r + Q*propBody;
        set(propHandles(i),'XData',propWorld(1,:), ...
                           'YData',propWorld(2,:), ...
                           'ZData',propWorld(3,:));
    end

    addpoints(traj,r(1),r(2),r(3));
    set(refMarker,'XData',rref(1),'YData',rref(2),'ZData',rref(3));
    set(statusText,'Position',[r(1),r(2),r(3)+0.25], ...
        'String',sprintf('v_{ref} = [%+.1f %+.1f %+.1f] m/s   yawdot_{ref} = %+.2f rad/s', ...
                         vrefN(1),vrefN(2),vrefN(3),yawRateRef));

    % Follow the vehicle with a fixed-size viewing box.
    xlim(ax,r(1) + [-followRadiusXY followRadiusXY]);
    ylim(ax,r(2) + [-followRadiusXY followRadiusXY]);
    zlo = max(-0.1,r(3)-followBelowZ);
    zhi = max(zlo+1.0,r(3)+followAboveZ);
    zlim(ax,[zlo zhi]);

    drawnow;

    % Pace the fixed-step simulation to wall-clock time.
    remaining = p.h - toc(ticStep);
    if remaining > 0
        pause(remaining);
    end
end

if ishandle(fig)
    close(fig);
end


%% Local helpers
function key_press(src,event)
    keys = getappdata(src,'keys');
    k = matlab.lang.makeValidName(lower(event.Key));
    keys.(k) = true;
    setappdata(src,'keys',keys);

    if strcmpi(event.Key,'escape') || strcmpi(event.Key,'q')
        setappdata(src,'quit',true);
    end
end

function key_release(src,event)
    keys = getappdata(src,'keys');
    k = matlab.lang.makeValidName(lower(event.Key));
    keys.(k) = false;
    setappdata(src,'keys',keys);
end

function tf = key_down(keys,key)
    k = matlab.lang.makeValidName(lower(key));
    tf = isfield(keys,k) && keys.(k);
end

function a = wrap_angle(a)
    a = atan2(sin(a),cos(a));
end

function A = numerical_jacobian(fun,x)
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
