function xdot = Qode(t,x,J)

Q = reshape(x(1:9),3,3);
w = x(10:12);

wdot = -J\(cross(w,J*w));
Qdot = Q*hat(w);

xdot = [Qdot(:); wdot];

end