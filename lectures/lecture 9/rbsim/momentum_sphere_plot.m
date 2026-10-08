%% momentum_sphere_plot.m
% Free rigid body in body angular momentum space.
%
% Integrates Euler's equations  dh/dt = h x (J \ h)  for a set of initial
% conditions and plots them on the momentum sphere together with the
% separatrices and equilibria.
%
% Requires on the MATLAB path:
%   plotMomentumSphere.m, plotRigidBodyTrajectory.m

clear; close all; clc;

%% ---------------- Parameters ----------------
Jprincipal = [1 2 3];      % principal moments (distinct)
rotateJ    = true;         % true: non-diagonal J (principal axes rotated)
hmag       = 1;            % angular momentum magnitude |h|
nPer       = 4;            % trajectories per family (each side of separatrix)
gap        = 0.06;         % angular offset from separatrix / equilibria [rad]
tFinal     = 80;           % integration time
opts       = odeset('RelTol',1e-10,'AbsTol',1e-12);

%% ---------------- Inertia matrix ----------------
J = diag(Jprincipal);
if rotateJ
    w = [0.4; -0.2; 0.3];                        % rotation vector
    W = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    R = expm(W);                                 % rotation matrix
    J = R*J*R.';
end

rhs = @(t,h) cross(h, J\h);                      % Euler's equations

%% ---------------- Momentum sphere ----------------
fig = figure('Color','w','Position',[100 100 950 820]);
ax  = axes(fig);
plotMomentumSphere(J, hmag, 'Parent', ax);

% Principal axes and moments (ascending), as used by the plot
V  = ax.UserData.rigidBody.V;
Jp = ax.UserData.rigidBody.Jp;

%% ---------------- Initial conditions ----------------
% Points on the great circle through v1 and v3, on both sides of the
% separatrix angle aSep, where tan^2(aSep) = (1/J1-1/J2)/(1/J2-1/J3)
aSep   = atan(sqrt((1/Jp(1)-1/Jp(2)) / (1/Jp(2)-1/Jp(3))));
alphas = [linspace(gap, aSep-gap, nPer), linspace(aSep+gap, pi/2-gap, nPer)];

H0 = zeros(3,0);
for a = alphas
    for s = [-1 1]
        if a < aSep    % orbits around +/- v1
            H0(:,end+1) = hmag*(s*cos(a)*V(:,1) + sin(a)*V(:,3)); %#ok<SAGROW>
        else           % orbits around +/- v3
            H0(:,end+1) = hmag*(cos(a)*V(:,1) + s*sin(a)*V(:,3)); %#ok<SAGROW>
        end
    end
end

%% ---------------- Integrate and plot trajectories ----------------
nIC  = size(H0,2);
traj = struct('t',cell(1,nIC),'h',cell(1,nIC));
for k = 1:nIC
    [t,h] = ode45(rhs, [0 tFinal], H0(:,k), opts);
    traj(k).t = t;  traj(k).h = h;

    plotRigidBodyTrajectory(ax, h, J);
end

rotate3d(fig, 'on');
