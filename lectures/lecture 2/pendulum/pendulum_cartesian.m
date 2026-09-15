function xdot = pendulum_cartesian(t,x)

g = 9.8;
m = 1.0;
length = 1.0;

r = x(1:2);
rdot = x(3:4);

T = m*(rdot'*rdot)/length - m*g*r(2);

rddot = [-T*r(1)/(m*length); (-T*r(2)/(m*length) - g)];

xdot = [rdot; rddot];

end