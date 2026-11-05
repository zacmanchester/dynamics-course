function dx = flyingcar_dynamics(x,u,p)
%FLYINGCAR_DYNAMICS Continuous-time flying-car rigid-body dynamics.
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% Input:
%   u = [T1; T2; T3; T4; T5; S; f]
%   four vertical cold-gas thrusters, the forward thruster, the suction
%   force and the driving force per wheel [N]
%
% Frames: inertial z up; body x forward, y left, z up.

    q = x(4:7);
    q = q/norm(q);
    v = x(8:10);
    w = x(11:13);

    Q = FlyingCarMath.qtoQ(q);

    rdot = Q*v;
    qdot = 0.5*FlyingCarMath.Lquat(q)*FlyingCarMath.Hmat()*w;

    T1 = u(1); T2 = u(2); T3 = u(3); T4 = u(4);
    T5 = u(5); S  = u(6); f  = u(7);

    % Tire forces (normal + driving) and their torque.
    [Ftire,tauTire] = tire_forces(x,f,p);

    % Total force: thrusters and suction + gravity + tires.
    F = [T5; 0; T1 + T2 + T3 + T4 - S] + Q'*[0; 0; -p.m*p.g] + Ftire;

    % Total torque: thrusters (lever arms b and a) + tires.
    tau = [p.b*(T1 + T4 - T2 - T3); ...
           p.a*(T3 + T4 - T1 - T2); ...
           0] + tauTire;

    vdot = F/p.m - FlyingCarMath.hat(w)*v;
    wdot = p.J \ (tau - FlyingCarMath.hat(w)*p.J*w);

    dx = [rdot; qdot; vdot; wdot];
end
