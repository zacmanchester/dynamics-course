function CL = lift_coefficient(alpha,p)
%LIFT_COEFFICIENT Static lift coefficient as a function of the angle of attack.
%   Linear below the stall angle, blending into the flat-plate value
%   2*sin^2(alpha)*cos(alpha) beyond it (stall). Works on arrays of alpha.

    sigma = 1./(1 + exp(-p.Ms*(alpha - p.alpha_s))) ...
          + 1./(1 + exp(-p.Ms*(-alpha - p.alpha_s)));
    CL = (1 - sigma).*(p.CL0 + p.CLa*alpha) ...
       + sigma.*2.*sign(alpha).*sin(alpha).^2.*cos(alpha);
end
