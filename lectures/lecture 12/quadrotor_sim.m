% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% where q = [scalar; vector] and maps body coordinates to world coordinates.

clear; clc; close all;
rng('default');

%% Quadrotor parameters
p.m  = 0.5;
p.ell = 0.1750;
p.J  = diag([0.0023, 0.0023, 0.004]);
p.g  = 9.81;
p.kt = 1.0;
p.km = 0.0245;
p.h  = 0.05;       % 20 Hz

Nx  = 13;          % quaternion state dimension
Nxr = 12;          % reduced/linearized state dimension
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
E0 = E(q0);
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

% Match the notebook's perturbed initial condition:
%   r = r0 + randn(3)
%   attitude Rodrigues parameter phi = [1;0;0]
xhist(:,1) = [r0 + randn(3,1); ...
              Lquat(q0)*rptoq([1;0;0]); ...
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

%% Optional simple animation
animate_quadrotor(xhist,thist,p.ell);


%% Local functions
function dx = quad_dynamics(x,u,p)
    r = x(1:3); %#ok<NASGU>
    q = x(4:7);
    q = q/norm(q);  % normalize just to be careful
    v = x(8:10);
    w = x(11:13);

    Q = qtoQ(q);

    rdot = Q*v;
    qdot = 0.5*Lquat(q)*Hmat()*w;

    thrustMap = [zeros(2,4); p.kt*ones(1,4)];
    vdot = Q'*[0;0;-p.g] + (1/p.m)*thrustMap*u - hat(w)*v;

    torqueMap = [0,          p.ell*p.kt, 0,             -p.ell*p.kt; ...
                -p.ell*p.kt, 0,          p.ell*p.kt,    0; ...
                 p.km,      -p.km,       p.km,          -p.km];

    wdot = p.J \ (-hat(w)*p.J*w + torqueMap*u);

    dx = [rdot; qdot; vdot; wdot];
end

function xn = quad_dynamics_rk4(x,u,p)
    % RK4 integration with zero-order hold on u.
    f1 = quad_dynamics(x,u,p);
    f2 = quad_dynamics(x + 0.5*p.h*f1,u,p);
    f3 = quad_dynamics(x + 0.5*p.h*f2,u,p);
    f4 = quad_dynamics(x + p.h*f3,u,p);

    xn = x + (p.h/6)*(f1 + 2*f2 + 2*f3 + f4);
    xn(4:7) = xn(4:7)/norm(xn(4:7));
end

function u = feedback_controller(x,x0,r0,v0,w0,uhover,K)
    q0 = x0(4:7);
    q  = x(4:7);

    % Relative attitude q_rel = q0^{-1} (*) q. For a unit quaternion,
    % L(q0)' is the left-multiplication matrix for q0^{-1}.
    phi = qtorp(Lquat(q0)'*q);

    dxr = [x(1:3)-r0; ...
           phi; ...
           x(8:10)-v0; ...
           x(11:13)-w0];

    u = uhover - K*dxr;
end

function M = E(q)
    % Map reduced 12-state perturbations into the 13-state quaternion space.
    M = blkdiag(eye(3),Gquat(q),eye(6));
end

function A = numerical_jacobian(fun,x)
    % Central-difference Jacobian.
    % Step size is scaled componentwise to balance truncation and roundoff.
    fx = fun(x);
    m = numel(fx);
    n = numel(x);
    A = zeros(m,n);

    baseStep = eps^(1/3);
    for i = 1:n
        hi = baseStep*max(1,abs(x(i)));
        xp = x; xm = x;
        xp(i) = xp(i) + hi;
        xm(i) = xm(i) - hi;
        A(:,i) = (fun(xp)-fun(xm))/(2*hi);
    end
end

function S = hat(v)
    S = [ 0,    -v(3),  v(2); ...
          v(3),  0,    -v(1); ...
         -v(2),  v(1),  0];
end

function M = Lquat(q)
    s = q(1);
    v = q(2:4);
    M = [s, -v'; ...
         v, s*eye(3) + hat(v)];
end

function H = Hmat()
    H = [zeros(1,3); eye(3)];
end

function Q = qtoQ(q)
    T = diag([1,-1,-1,-1]);
    H = Hmat();
    Q = H'*T*Lquat(q)*T*Lquat(q)*H;
end

function G = Gquat(q)
    G = Lquat(q)*Hmat();
end

function q = rptoq(phi)
    q = [1; phi]/sqrt(1 + phi'*phi);
end

function phi = qtorp(q)
    phi = q(2:4)/q(1);
end

function animate_quadrotor(xhist,thist,armLength)
    fig = figure('Name','Quadrotor animation');
    ax = axes(fig);
    hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    xlabel(ax,'x [m]'); ylabel(ax,'y [m]'); zlabel(ax,'z [m]');
    title(ax,'Quadrotor animation');

    r = xhist(1:3,:);
    pad = max(0.5,2*armLength);
    xlim(ax,[min(r(1,:))-pad, max(r(1,:))+pad]);
    ylim(ax,[min(r(2,:))-pad, max(r(2,:))+pad]);
    zlim(ax,[min(0,min(r(3,:))-pad), max(r(3,:))+pad]);
    view(ax,3);

    plot3(ax,r(1,:),r(2,:),r(3,:),'--');
    arm1 = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'-','LineWidth',2);
    arm2 = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'-','LineWidth',2);

    bodyPts1 = armLength*[1, -1; 0, 0; 0, 0];
    bodyPts2 = armLength*[0, 0; 1, -1; 0, 0];

    for k = 1:size(xhist,2)
        rk = xhist(1:3,k);
        qk = xhist(4:7,k);
        qk = qk/norm(qk);
        Qk = qtoQ(qk);

        p1 = rk + Qk*bodyPts1;
        p2 = rk + Qk*bodyPts2;

        set(arm1,'XData',p1(1,:),'YData',p1(2,:),'ZData',p1(3,:));
        set(arm2,'XData',p2(1,:),'YData',p2(2,:),'ZData',p2(3,:));
        drawnow;

        if k < size(xhist,2)
            pause(max(0,thist(k+1)-thist(k)));
        end
    end
end
