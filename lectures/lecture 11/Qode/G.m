function Gmat = G(q)
    s = q(1);
    v = q(2:4);
    Gmat = [-v'; s*eye(3)+hat(v)];
end