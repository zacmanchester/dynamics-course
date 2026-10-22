function xn = quad_dynamics_rk4(x,u,p)
%QUAD_DYNAMICS_RK4 One RK4 step with zero-order hold on the input.

    f1 = quad_dynamics(x,u,p);
    f2 = quad_dynamics(x + 0.5*p.h*f1,u,p);
    f3 = quad_dynamics(x + 0.5*p.h*f2,u,p);
    f4 = quad_dynamics(x + p.h*f3,u,p);

    xn = x + (p.h/6)*(f1 + 2*f2 + 2*f3 + f4);
    xn(4:7) = xn(4:7)/norm(xn(4:7));
end
