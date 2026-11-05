function p = flyingcar_params()
%FLYINGCAR_PARAMS Parameters of the flying car, shared by all demos.
%
% Frames: inertial z up; body x forward, y left, z up.
% Thrusters T1..T4 point along body z and sit at the four wheels:
%   T1 front-left (+a,+b), T2 front-right (+a,-b),
%   T3 rear-right (-a,-b),  T4 rear-left  (-a,+b).
% T5 points forward along body x. The suction fans pull along -z.

    p.m   = 1500;                       % mass [kg]
    p.J   = diag([600, 2500, 2800]);    % inertia [kg m^2]
    p.g   = 9.81;
    p.h   = 0.01;                       % 100 Hz

    p.a   = 1.4;                        % half wheelbase [m]
    p.b   = 0.9;                        % half track [m]
    p.hw  = 0.45;                       % center of mass height above the wheels [m]

    % Wheel (and thruster) positions in the body frame.
    p.rw = [ p.a,  p.a, -p.a, -p.a; ...
             p.b, -p.b, -p.b,  p.b; ...
            -p.hw, -p.hw, -p.hw, -p.hw];

    % Actuator limits.
    p.Tmax  = p.m*p.g;                  % each vertical thruster [N]
    p.T5max = p.m*p.g;                  % forward thruster [N]
    p.Smax  = 1.5*p.m*p.g;              % suction [N]

    % Tires: a spring-damper against the ground, and friction-limited
    % driving force.
    p.kg = 2e5;                         % tire stiffness [N/m]
    p.cg = 1e4;                         % tire damping [N s/m]
    p.mu = 1.0;                         % friction coefficient
    p.cs = 2e4;                         % sideways friction per unit slip speed [N s/m]

    % Flat ground by default; a demo may switch the loop on.
    p.loop.on = false;
end
