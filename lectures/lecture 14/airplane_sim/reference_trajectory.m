function [xref,tphase] = reference_trajectory(t,p)
%REFERENCE_TRAJECTORY Racetrack reference state at time t.
%
%   xref   = [r_N(3); q(4); v_B(3); omega_B(3)]
%   tphase = end time of each phase
%
% Left-hand racetrack. Take-off and landing on the same runway along +x:
%   1 ground roll   2 climb   3 left turn   4 downwind leg
%   5 left turn     6 final approach        7 braking

    V  = p.Vc;                            % cruise speed
    a  = p.a_runway;                      % take-off acceleration
    ab = p.a_brake;                       % braking deceleration
    R  = V^2/(p.g*tan(p.bank));           % turn radius for the chosen bank angle
    psidot = V/R;                         % turn rate

    % Phase durations and end times.
    t1 = V/a;                             % ground roll
    t2 = p.hc/(V*sin(p.gamma_up));        % climb
    t3 = pi/psidot;                       % 180 deg turn
    t4 = p.downwind/V;                    % downwind leg
    t5 = t3;                              % 180 deg turn
    t6 = p.hc/(V*sin(p.gamma_down));      % final approach
    t7 = V/ab;                            % braking
    tphase = cumsum([t1 t2 t3 t4 t5 t6 t7]);

    % Key points along the runway axis.
    x1 = 0.5*a*t1^2;                      % lift-off
    x2 = x1 + V*cos(p.gamma_up)*t2;       % top of climb, start of turn 1
    x4 = x2 - p.downwind;                 % end of downwind leg, start of turn 2
    x6 = x4 + V*cos(p.gamma_down)*t6;     % touchdown

    % Defaults: cruise speed, level path, heading +x, no bank, no turning.
    speed = V; gamma = 0; psi = 0; bank = 0; psid = 0; onGround = false;

    if t < tphase(1)                                % 1 ground roll
        onGround = true;
        speed = a*t;
        r = [0.5*a*t^2; 0; 0];
    elseif t < tphase(2)                            % 2 climb
        s = t - tphase(1);
        gamma = p.gamma_up;
        r = [x1 + V*cos(gamma)*s; 0; V*sin(gamma)*s];
    elseif t < tphase(3)                            % 3 left turn (centre at (x2,R))
        s = t - tphase(2);
        psi = psidot*s;  psid = psidot;  bank = -p.bank;
        r = [x2 + R*sin(psi); R - R*cos(psi); p.hc];
    elseif t < tphase(4)                            % 4 downwind leg, heading -x
        s = t - tphase(3);
        psi = pi;
        r = [x2 - V*s; 2*R; p.hc];
    elseif t < tphase(5)                            % 5 left turn (centre at (x4,R))
        s = t - tphase(4);
        psi = pi + psidot*s;  psid = psidot;  bank = -p.bank;
        r = [x4 + R*sin(psi); R - R*cos(psi); p.hc];
    elseif t < tphase(6)                            % 6 final approach
        s = t - tphase(5);
        gamma = -p.gamma_down;
        r = [x4 + V*cos(gamma)*s; 0; p.hc + V*sin(gamma)*s];
    elseif t < tphase(7)                            % 7 braking
        onGround = true;
        s = t - tphase(6);
        speed = V - ab*s;
        r = [x6 + V*s - 0.5*ab*s^2; 0; 0];
    else                                            % stopped on the runway
        onGround = true;
        speed = 0;
        r = [x6 + 0.5*V*t7; 0; 0];
    end

    % Reference attitude: yaw psi, then pitch (trim pitch + path angle), then
    % roll. On the ground the airplane is level. In this frame a nose-up pitch
    % is a negative rotation about body y, and a left turn needs the left wing
    % down, i.e. a negative roll.
    if onGround
        theta = 0;
    else
        theta = p.theta0 + gamma;
    end
    qyaw   = [cos(psi/2);   0; 0; sin(psi/2)];
    qpitch = [cos(theta/2); 0; -sin(theta/2); 0];
    qroll  = [cos(bank/2);  sin(bank/2); 0; 0];
    qref = AirplaneMath.Lquat(qyaw)*AirplaneMath.Lquat(qpitch)*qroll;
    Qref = AirplaneMath.qtoQ(qref);

    % Reference velocity along the path and angular velocity of the turn,
    % both expressed in the body frame.
    vN   = speed*[cos(psi)*cos(gamma); sin(psi)*cos(gamma); sin(gamma)];
    vref = Qref'*vN;
    wref = Qref'*[0; 0; psid];

    xref = [r; qref; vref; wref];
end
