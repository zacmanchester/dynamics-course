function xdot = quatode(t,x,J)
    q = x(1:4);
    w = x(5:7);

    kp = 1;

    wdot = -J\(cross(w,J*w));
    qdot = G(q)*w;
    %qdot = G(q)*w - kp*(q'*q-1)*q;

    xdot = [qdot; wdot];
end