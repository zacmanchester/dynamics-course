function xdot = quatode(t,x,J)
    q = x(1:4);
    q = q/norm(q); %renormalize
    w = x(5:7);

    wdot = -J\(cross(w,J*w));
    qdot = G(q)*w;

    xdot = [qdot; wdot];
end
