function xn = rocket_dynamics_rk4(x,u,p)
%ROCKET_DYNAMICS_RK4 One RK4 step with zero-order hold on the input.

    f1 = rocket_dynamics(x,u,p);
    f2 = rocket_dynamics(x + 0.5*p.h*f1,u,p);
    f3 = rocket_dynamics(x + 0.5*p.h*f2,u,p);
    f4 = rocket_dynamics(x + p.h*f3,u,p);

    xn = x + (p.h/6)*(f1 + 2*f2 + 2*f3 + f4);
    xn(4:7) = xn(4:7)/norm(xn(4:7));

    % Landing pad: the legs touch the ground when the center of mass is at
    % z = ell. The pad holds the rocket there and the legs stop all motion.
    if xn(3) < p.ell
        xn(3) = p.ell;
        xn(8:13) = 0;
    end
end
