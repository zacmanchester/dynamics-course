function u = feedback_controller(x,xref,utrim,K,umin,umax)
%FEEDBACK_CONTROLLER Reduced-coordinate LQR tracking controller.
%   The state error is expressed relative to the reference state xref, so the
%   gain K computed about the straight-and-level trim can be reused along the
%   whole trajectory.

    rref = xref(1:3);
    qref = xref(4:7);
    vref = xref(8:10);
    wref = xref(11:13);
    Qref = AirplaneMath.qtoQ(qref);

    q = x(4:7)/norm(x(4:7));

    % Attitude error q_rel = qref^{-1} (*) q.
    phi = AirplaneMath.qtorp(AirplaneMath.Lquat(qref)'*q);

    dxr = [Qref'*(x(1:3)-rref); ...
           phi; ...
           x(8:10)-vref; ...
           x(11:13)-wref];

    u = utrim - K*dxr;

    % Actuator limits: thrust in [0,Tmax], surfaces within +/- dmax.
    u = min(max(u,umin),umax);
end
