function dx = rocket_dynamics(x,u,p)
%ROCKET_DYNAMICS Continuous-time thrust-vectored rocket rigid-body dynamics.
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% Input:
%   u = [T; alpha; beta]
%   thrust magnitude [N], gimbal angles about body y and body x [rad]
%
% Frames: inertial z up; body z along the rocket axis (nose up).
% The engine sits a distance ell below the center of mass.

    q = x(4:7);
    q = q/norm(q);
    v = x(8:10);
    w = x(11:13);

    Q = RocketMath.qtoQ(q);

    rdot = Q*v;
    qdot = 0.5*RocketMath.Lquat(q)*RocketMath.Hmat()*w;

    T     = u(1);
    alpha = u(2);
    beta  = u(3);

    % Thrust force in the body frame (gimballed about y by alpha, then x by beta).
    FT = T*[sin(alpha); sin(beta)*cos(alpha); cos(beta)*cos(alpha)];

    % Total force: thrust + gravity.
    F = FT + Q'*[0; 0; -p.m*p.g];

    % Moment of the thrust about the center of mass.
    rT  = [0; 0; -p.ell];
    tau = cross(rT,FT);

    vdot = F/p.m - RocketMath.hat(w)*v;
    wdot = p.J \ (tau - RocketMath.hat(w)*p.J*w);

    dx = [rdot; qdot; vdot; wdot];
end
