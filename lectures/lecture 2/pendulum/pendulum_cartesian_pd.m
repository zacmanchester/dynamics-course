function xdot = pendulum_cartesian_pd(t,x)

g = 9.8;
m = 1.0;
length = 1.0;

r = x(1:2);
rdot = x(3:4);

c = (sqrt(r'*r) - length);
cdot = r'*rdot/sqrt(r'*r);
kp = 10;
kd = 1;

T = m*(rdot'*rdot)/length - m*g*r(2) + kp*c + kd*cdot;

rddot = [-T*r(1)/(m*length); (-T*r(2)/(m*length) - g)];

xdot = [rdot; rddot];

end