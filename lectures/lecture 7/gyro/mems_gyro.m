function xdot = mems_gyro(t,x,omega,m,k,c,a)

[wz, wzdot] = omega(t);

S = hat([0; 0; wz]);
Sdot = hat([0; 0; wzdot]);

w0 = sqrt(k/m);

r = x(1:3);
rdot = x(4:6);

rddot = -(2*S*rdot + Sdot*r + S*S*r + c*rdot + (k/m)*r + (c/m)*rdot - [a*sin(w0*t); 0; 0]);

xdot = [rdot; rddot];

end