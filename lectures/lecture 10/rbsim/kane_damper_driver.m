%% kane_damper_driver.m
% Animate a box-shaped rigid body with a spherical (Kane) damper.
%
% Model: a rigid body containing a spherical cavity. A sphere with
% isotropic inertia j sits in the cavity at the center of mass, and viscous
% fluid couples its rotation to the body with damping coefficient c. The
% fluid dissipates energy while total angular momentum is conserved, so the
% body drifts from spin about its minor axis (maximum energy for a given
% |h|) to spin about its major axis (minimum energy).
%
% Body-frame equations (omega: body rate, omegaD: damper rate, both
% absolute rates expressed in body coordinates; J excludes the damper's
% rotational inertia j):
%   J  dw/dt  + w x (J w)  =  c (wD - w)
%   j  dwD/dt + j (w x wD) = -c (wD - w)
%   dq/dt = (1/2) q (x) [0; w]
%
% Total body-frame angular momentum h = J w + j wD has constant magnitude.
% The momentum sphere panel uses the locked inertia Jlock = J + j*eye(3),
% whose energy ellipsoids and separatrices describe the body once the
% damper has stopped slipping.
%
% Requires on the MATLAB path:  makeBoxBody.m, animateRigidBody.m,
%   plotMomentumSphere.m, plotRigidBodyTrajectory.m

clear; close all; clc;

%% ---------------- Body ----------------
body = makeBoxBody('Lx',1.2, 'Ly',0.8, 'Lz',0.4);   % distinct sides -> distinct inertias
J    = body.J;

[~, order] = sort(diag(J));
iMin = order(1);  iMax = order(3);
axisNames = 'xyz';
fprintf('Principal moments: [%.4g %.4g %.4g]  (minor axis: %c, major axis: %c)\n', ...
        diag(J), axisNames(iMin), axisNames(iMax));

%% ---------------- Damper ----------------
jD = 0.10*max(diag(J));    % damper sphere inertia
cD = 0.10*max(diag(J));    % viscous coupling coefficient [torque per rad/s]
Jlock = J + jD*eye(3);     % inertia with the damper locked to the body

%% ---------------- Initial conditions ----------------
Omega   = pi;            % nominal spin rate [rad/s] (1 rev/s)
perturb = 0.02;            % relative perturbation off the spin axis
spinAxis = iMin;           % start spinning about the minor axis

w0 = perturb*Omega*ones(3,1);
w0(spinAxis) = Omega;
wD0 = w0;                  % damper initially moving with the body
q0  = [1; 0; 0; 0];        % body frame aligned with inertial frame

%% ---------------- Integrate ----------------
tFinal  = 120;              % [s]
nFrames = 2000;            % output samples = animation frames
tspan   = linspace(0, tFinal, nFrames);
opts    = odeset('RelTol',1e-10, 'AbsTol',1e-12);

[t, Y] = ode45(@(t,y) kaneDamperODE(y, J, jD, cD), tspan, [w0; wD0; q0], opts);

%% ---------------- Post-process ----------------
N  = numel(t);
W  = Y(:,1:3).';                  % body rate (body frame), 3xN
WD = Y(:,4:6).';                  % damper rate (body frame), 3xN
hB = J*W + jD*WD;                 % total body-frame angular momentum
T  = 0.5*sum(W.*(J*W), 1) + 0.5*jD*sum(WD.^2, 1);   % total kinetic energy

R  = zeros(3,3,N);                % body -> inertial rotations
wI = zeros(3,N);                  % inertial body angular velocity
hI = zeros(3,N);                  % inertial angular momentum (constant)
for k = 1:N
    q = Y(k,7:10).';
    R(:,:,k) = quatToRotm(q/norm(q));
    wI(:,k)  = R(:,:,k)*W(:,k);
    hI(:,k)  = R(:,:,k)*hB(:,k);
end
fprintf('Max relative drift of inertial h: %.2e\n', ...
        max(vecnorm(hI - hI(:,1), 2, 1))/norm(hI(:,1)));
fprintf('Energy: T(0) = %.4g  ->  T(end) = %.4g  (min possible %.4g)\n', ...
        T(1), T(end), norm(hB(:,1))^2/(2*max(diag(Jlock))));

%% ---------------- Animate ----------------
% Top left: box in the inertial frame.  Top right: body-frame h on the
% momentum sphere (locked inertia), trace colored by the decaying energy.
% Bottom: energy vs time, with the min / separatrix / max energy levels.
animateRigidBody(t, R, body, 'h', hI(:,1), 'omega', wI, ...
                 'hBody', hB, 'J', Jlock, 'Energy', T, ...
                 'ShowEnergy', true, 'Speed', 3);
% To save a video instead, add:  'VideoFile', 'kane_damper.mp4'

%% ======================= Local functions =======================
function dy = kaneDamperODE(y, J, jD, cD)
% Rigid body with a spherical (Kane) damper, plus quaternion kinematics.
w  = y(1:3);
wD = y(4:6);
q  = y(7:10);
tauD = cD*(wD - w);                          % viscous torque on the body
dw   = J \ (tauD - cross(w, J*w));
dwD  = -tauD/jD - cross(w, wD);
q0 = q(1);  qv = q(2:4);
dq = 0.5*[-qv.'*w; q0*w + cross(qv, w)];     % (1/2) q (x) [0; w]
dy = [dw; dwD; dq];
end

function R = quatToRotm(q)
% Rotation matrix (body -> inertial) from a unit quaternion [q0; q1; q2; q3].
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
R = [1-2*(q2^2+q3^2),   2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
       2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2),   2*(q2*q3-q0*q1);
       2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
end
