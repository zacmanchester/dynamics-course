function dx = airplane_dynamics(x,u,p)
%AIRPLANE_DYNAMICS Continuous-time fixed-wing airplane rigid-body dynamics.
%
% State:
%   x = [r_N(3); q(4); v_B(3); omega_B(3)]
% Input:
%   u = [T; delta_a; delta_e; delta_r]
%   thrust [N], aileron, elevator and rudder deflections [rad]
%
% Frames: inertial z up; body x forward, y toward the left wing, z up.
% The Newton-Euler equations are the same as for the quadrotor. Only the
% external force F and moment tau change.

    q = x(4:7);
    q = q/norm(q);
    v = x(8:10);
    w = x(11:13);

    Q = AirplaneMath.qtoQ(q);

    rdot = Q*v;
    qdot = 0.5*AirplaneMath.Lquat(q)*AirplaneMath.Hmat()*w;

    T  = u(1);
    da = u(2);
    de = u(3);
    dr = u(4);

    % Air-relative velocity (no wind): airspeed, angle of attack, sideslip.
    Va    = norm(v);
    alpha = atan2(-v(3), v(1));
    beta  = atan2(v(2), sqrt(v(1)^2 + v(3)^2));
    qbar  = 0.5*p.rho*Va^2;

    % Non-dimensional angular rates (Vs avoids dividing by zero at rest).
    Vs   = max(Va,1);
    phat = p.b/(2*Vs)*w(1);
    qhat = p.c/(2*Vs)*w(2);
    rhat = p.b/(2*Vs)*w(3);

    % Aerodynamic coefficients (stability-derivative approximation).
    % The static lift curve CL(alpha), including stall, is in lift_coefficient.
    CL = lift_coefficient(alpha,p) + p.CLq*qhat + p.CLde*de;
    CD = p.CD0 + CL^2/(pi*p.e*p.AR);
    CY = p.CYb*beta + p.CYda*da + p.CYdr*dr;
    Cl = p.Clb*beta + p.Clp*phat + p.Clr*rhat + p.Clda*da + p.Cldr*dr;
    Cm = p.Cm0 + p.Cma*alpha + p.Cmq*qhat + p.Cmde*de;
    Cn = p.Cnb*beta + p.Cnp*phat + p.Cnr*rhat + p.Cnda*da + p.Cndr*dr;

    % Aerodynamic forces.
    L = qbar*p.S*CL;
    D = qbar*p.S*CD;
    Y = qbar*p.S*CY;

    % Total force in the body frame: thrust + gravity + aerodynamics.
    F = [T; 0; 0] ...
      + Q'*[0; 0; -p.m*p.g] ...
      + [L*sin(alpha) - D*cos(alpha); Y; L*cos(alpha) + D*sin(alpha)];

    % Aerodynamic moments about the center of mass.
    tau = qbar*p.S*[p.b*Cl; p.c*Cm; p.b*Cn];

    vdot = F/p.m - AirplaneMath.hat(w)*v;
    wdot = p.J \ (tau - AirplaneMath.hat(w)*p.J*w);

    dx = [rdot; qdot; vdot; wdot];
end
