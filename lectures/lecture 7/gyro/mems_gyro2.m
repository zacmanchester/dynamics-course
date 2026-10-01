function xdot = mems_gyro2(t,x,omega,m,k,c,a)

[wz, wzdot] = omega(t);

w0 = sqrt(k/m);

rx = x(1);
ry = x(2);
rxdot = x(3);
rydot = x(4);

rxddot = -(k/m)*rx - (c/m)*rxdot + a*sin(w0*t);

ryddot = -(k/m)*ry - (c/m)*rydot - 2*wz*rxdot;


xdot = [rxdot; rydot; rxddot; ryddot];

end