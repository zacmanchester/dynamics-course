function xdot = ant_record(t,x,omega,c)

S = hat(omega);

r = x(1:3);
rdot = x(4:6);

rddot = -(2*S*rdot + S*S*r + c*rdot);

xdot = [rdot; rddot];

end