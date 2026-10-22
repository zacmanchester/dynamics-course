function dx = quad_dynamics(x,u,p)
%QUAD_DYNAMICS Continuous-time quadrotor rigid-body dynamics.
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]

    q = x(4:7);
    q = q/norm(q);
    v = x(8:10);
    w = x(11:13);

    Q = QuadrotorMath.qtoQ(q);

    rdot = Q*v;
    qdot = 0.5*QuadrotorMath.Lquat(q)*QuadrotorMath.Hmat()*w;

    thrustMap = [zeros(2,4); p.kt*ones(1,4)];
    vdot = Q'*[0;0;-p.g] + (1/p.m)*thrustMap*u ...
           - QuadrotorMath.hat(w)*v;

    torqueMap = [0,             p.ell*p.kt, 0,             -p.ell*p.kt; ...
                -p.ell*p.kt,   0,          p.ell*p.kt,     0; ...
                 p.km,         -p.km,      p.km,           -p.km];

    wdot = p.J \ (-QuadrotorMath.hat(w)*p.J*w + torqueMap*u);

    dx = [rdot; qdot; vdot; wdot];
end
