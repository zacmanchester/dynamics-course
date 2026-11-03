function xn = airplane_dynamics_rk4(x,u,p)
%AIRPLANE_DYNAMICS_RK4 One RK4 step with zero-order hold on the input.

    f1 = airplane_dynamics(x,u,p);
    f2 = airplane_dynamics(x + 0.5*p.h*f1,u,p);
    f3 = airplane_dynamics(x + 0.5*p.h*f2,u,p);
    f4 = airplane_dynamics(x + p.h*f3,u,p);

    xn = x + (p.h/6)*(f1 + 2*f2 + 2*f3 + f4);
    xn(4:7) = xn(4:7)/norm(xn(4:7));

    % Runway: the ground holds the airplane at z = 0, removes any downward
    % velocity and slows it down through friction. The landing gear keeps the
    % wings level and the nose at or above the horizon (the airplane can still
    % rotate nose-up for take-off), and prevents any roll or yaw rate.
    if xn(3) < 0
        xn(3) = 0;

        Q  = AirplaneMath.qtoQ(xn(4:7));
        vN = Q*xn(8:10);
        vN(3) = max(vN(3),0);

        psi   = atan2(Q(2,1),Q(1,1));                  % heading
        theta = max(asin(Q(3,1)),0);                   % pitch, nose up >= 0
        qyaw   = [cos(psi/2);   0; 0; sin(psi/2)];
        qpitch = [cos(theta/2); 0; -sin(theta/2); 0];
        xn(4:7) = AirplaneMath.Lquat(qyaw)*qpitch;
        Q = AirplaneMath.qtoQ(xn(4:7));

        xn(8:10) = Q'*vN;
        xn(8) = max(xn(8) - p.h*p.mu*p.g, 0);

        xn(11) = 0;
        xn(13) = 0;
        if theta == 0
            xn(12) = min(xn(12),0);
        end
    end
end
