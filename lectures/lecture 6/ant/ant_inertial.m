function xdot = ant_inertial(t,x,omega,c)

S = hat(omega);

r = x(1:3);
rdot = x(4:6);

rddot = [0; 0; 0] - c*(rdot - S*r);

xdot = [rdot; rddot];

end