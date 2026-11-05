function u = feedback_controller(x,xref,uref,K,umin,umax)
%FEEDBACK_CONTROLLER Reduced-coordinate LQR tracking controller.
%   The state error is expressed relative to the reference state xref and
%   the feedback is added to the feedforward input uref.

    rref = xref(1:3);
    qref = xref(4:7);
    vref = xref(8:10);
    wref = xref(11:13);
    Qref = FlyingCarMath.qtoQ(qref);

    q = x(4:7)/norm(x(4:7));

    % Attitude error q_rel = qref^{-1} (*) q.
    phi = FlyingCarMath.qtorp(FlyingCarMath.Lquat(qref)'*q);

    dxr = [Qref'*(x(1:3)-rref); ...
           phi; ...
           x(8:10)-vref; ...
           x(11:13)-wref];

    u = uref - K*dxr;

    % Actuator limits.
    u = min(max(u,umin),umax);
end
