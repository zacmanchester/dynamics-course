%% tumbling_T_driver.m
% Animate a tumbling T-shaped rigid body (Dzhanibekov / tennis racket effect).
%
% State integrated with ode45:
%   y = [h; q]
%     h : body-frame angular momentum (3x1),  dh/dt = h x omega,  omega = J \ h
%     q : unit quaternion [q0; qv] (scalar first), body -> inertial,
%         dq/dt = (1/2) q (x) [0; omega]
%
% Requires on the MATLAB path:  makeTBody.m, animateRigidBody.m,
%   plotMomentumSphere.m, plotRigidBodyTrajectory.m (momentum sphere panel)

clear; close all; clc;

%% ---------------- Body ----------------
body = makeTBody('BarLength',1.3, 'StemLength',0.8, 'Width',0.15, 'Thickness',0.15);
J    = body.J;

% Principal moments are on the diagonal of J (body axes are principal axes
% for the T). Identify the intermediate axis.
[~, order] = sort(diag(J));
iMin = order(1);
iMid = order(2);
iMaj = order(3);
axisNames = 'xyz';
fprintf('Principal moments: [%.4g %.4g %.4g]  (intermediate axis: %c)\n', ...
        diag(J), axisNames(iMid));

%% ---------------- Initial conditions ----------------
Omega   = 2*pi;        % nominal spin rate [rad/s] (1 rev/s)
perturb = 0.3;        % relative perturbation off the spin axis
spinAxis = iMaj;       % spin about the intermediate (unstable) axis

w0 = perturb*Omega*ones(3,1);
w0(spinAxis) = Omega;
h0 = J*w0;             % body angular momentum
q0 = [1; 0; 0; 0];     % body frame initially aligned with inertial frame

%% ---------------- Integrate ----------------
tFinal  = 30;          % [s]
nFrames = 900;         % output samples (uniform in time) = animation frames
tspan   = linspace(0, tFinal, nFrames);
opts    = odeset('RelTol',1e-10, 'AbsTol',1e-12);

[t, Y] = ode45(@(t,y) rigidBodyODE(y, J), tspan, [h0; q0], opts);

%% ---------------- Post-process ----------------
N  = numel(t);
hB = Y(:,1:3).';                   % body-frame angular momentum, 3xN
R  = zeros(3,3,N);                 % body -> inertial rotations
wI = zeros(3,N);                   % inertial angular velocity
hI = zeros(3,N);                   % inertial angular momentum (should be constant)
for k = 1:N
    q = Y(k,4:7).';
    R(:,:,k) = quatToRotm(q/norm(q));
    wI(:,k)  = R(:,:,k)*(J\hB(:,k));
    hI(:,k)  = R(:,:,k)*hB(:,k);
end
fprintf('Max relative drift of inertial h: %.2e\n', ...
        max(vecnorm(hI - hI(:,1), 2, 1))/norm(hI(:,1)));

%% ---------------- Animate ----------------
% Left: T in the inertial frame.  Right: body-frame h on the momentum sphere.
animateRigidBody(t, R, body, 'h', hI(:,1), 'omega', wI, ...
                 'hBody', hB, 'J', J, 'Speed', 1);
% To save a video instead, add:  'VideoFile', 'tumbling_T.mp4'

%% ======================= Local functions =======================
function dy = rigidBodyODE(y, J)
% Euler's equations plus quaternion kinematics.
h  = y(1:3);
q  = y(4:7);
w  = J\h;                                   % body angular velocity
dh = cross(h, w);
q0 = q(1);  qv = q(2:4);
dq = 0.5*[-qv.'*w; q0*w + cross(qv, w)];    % (1/2) q (x) [0; w]
dy = [dh; dq];
end

function R = quatToRotm(q)
% Rotation matrix (body -> inertial) from a unit quaternion [q0; q1; q2; q3].
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
R = [1-2*(q2^2+q3^2),   2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
       2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2),   2*(q2*q3-q0*q1);
       2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
end
