function u = feedback_controller(x,x0,r0,v0,w0,uhover,K)
%FEEDBACK_CONTROLLER Reduced-coordinate LQR feedback controller.

    q0 = x0(4:7);
    q  = x(4:7);

    % Attitude error q_rel = q0^{-1} (*) q.
    phi = QuadrotorMath.qtorp(QuadrotorMath.Lquat(q0)'*q);

    dxr = [x(1:3)-r0; ...
           phi; ...
           x(8:10)-v0; ...
           x(11:13)-w0];

    u = uhover - K*dxr;
end
