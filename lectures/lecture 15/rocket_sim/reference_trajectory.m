function xref = reference_trajectory(t,p)
%REFERENCE_TRAJECTORY Landing reference state at time t.
%
%   xref = [r_N(3); q(4); v_B(3); omega_B(3)]
%
% Upright, above the pad, descending from h0 to the standing height ell along
% a smooth profile that starts and ends with zero vertical speed.

    s = min(t/p.Tland,1);
    z  = p.ell + (p.h0 - p.ell)*(1 - 3*s^2 + 2*s^3);
    vz = (p.h0 - p.ell)*(-6*s + 6*s^2)/p.Tland;

    xref = [0; 0; z; ...          % position: above the pad
            1; 0; 0; 0; ...       % attitude: upright
            0; 0; vz; ...         % body velocity (body z is up when upright)
            0; 0; 0];             % angular velocity
end
