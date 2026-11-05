function [Ftire,tauTire,N] = tire_forces(x,f,p)
%TIRE_FORCES Tire forces and torque in the body frame.
%   Each wheel i pushes with a normal force N_i along the ground normal n_i
%   and a driving force f_i along the forward tangent t_i (world frame):
%       F_i = N_i n_i + f_i t_i
%   The body-frame force is Q'F_i and its torque is r_i x Q'F_i.
%   N_i comes from a spring-damper tire, f_i is the commanded driving force
%   f limited by friction, and both are zero when the wheel is in the air.
%   On top of the notes, a sideways friction force stops the tire from
%   sliding sideways (otherwise the car would skate after landing).

    r = x(1:3);
    q = x(4:7)/norm(x(4:7));
    v = x(8:10);
    w = x(11:13);
    Q = FlyingCarMath.qtoQ(q);

    Ftire   = zeros(3,1);
    tauTire = zeros(3,1);
    N       = zeros(4,1);

    for i = 1:4
        ri = p.rw(:,i);                     % wheel position, body frame
        pw = r + Q*ri;                      % wheel position, world frame
        vw = Q*(v + FlyingCarMath.hat(w)*ri);   % wheel velocity, world frame

        [pen,n] = ground(pw,p);
        if pen <= 0
            continue;                       % wheel in the air
        end

        % Normal force: spring + damper, never pulling.
        N(i) = max(p.kg*pen - p.cg*(n'*vw), 0);

        % Driving force along the forward direction, limited by friction.
        t  = Q*[1; 0; 0];
        t  = t - (n'*t)*n;
        t  = t/norm(t);
        fi = min(max(f, -p.mu*N(i)), p.mu*N(i));

        % Sideways friction: opposes the sideways velocity of the wheel.
        s  = cross(n,t);
        fs = -min(max(p.cs*(s'*vw), -p.mu*N(i)), p.mu*N(i));

        Fi = Q'*(N(i)*n + fi*t + fs*s);
        Ftire   = Ftire + Fi;
        tauTire = tauTire + cross(ri,Fi);
    end
end
