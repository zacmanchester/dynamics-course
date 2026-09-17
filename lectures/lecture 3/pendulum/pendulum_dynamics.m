function xdot = pendulum_dynamics(t,x)

g = 9.8;
length = 1.0;

theta = x(1);
thetadot = x(2);

thetaddot = -(g/length)*sin(theta);

xdot = [thetadot; thetaddot];

end