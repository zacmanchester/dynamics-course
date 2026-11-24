%% Jacobian and eigenvalues of the linearized CR3BP at L1 and L4
mu = 0.012150585609624;               % Earth-Moon mass ratio

% --- L1 ---
% The Lagrange points are equilibria: places where a spacecraft at rest in
% the rotating frame stays at rest. With zero velocity, the equations of
% motion reduce to grad U = 0. The collinear points L1, L2, and L3 lie on
% the x axis (y = z = 0), where dU/dy = dU/dz = 0 automatically, so we only
% need to solve dU/dx = 0, a single equation in x.
%
% On the x axis, r1 = |x + mu| and r2 = |x - 1 + mu|, so
%     dU/dx = -x + (1-mu)(x+mu)/|x+mu|^3 + mu(x-1+mu)/|x-1+mu|^3.
% (The absolute values make the same formula valid on either side of each
% primary.)
dUdx_axis = @(x) -x + (1-mu)*(x+mu)./abs(x+mu).^3 + mu*(x-1+mu)./abs(x-1+mu).^3;

% L1 lies between the Earth (x = -mu) and the Moon (x = 1-mu). fzero needs
% an interval where dU/dx changes sign. dU/dx blows up near each primary,
% so we start the bracket a bit away from the Moon: on [0.5, 1-mu-0.01],
% dU/dx changes sign exactly once, at L1.
x_L1 = fzero(dUdx_axis, [0.5, 1 - mu - 0.01]);
L1 = [x_L1; 0; 0];

% --- L4 ---
% L4 forms an equilateral triangle with the two primaries, 60 degrees ahead
% of the Moon in its orbit.
L4 = [0.5 - mu; sqrt(3)/2; 0];

% --- Linearize at each point ---
Omega = [0 2 0; -2 0 0; 0 0 0];       % Coriolis terms
pts = {L1, L4};  names = {'L1', 'L4'};
for k = 1:numel(pts)
    r = pts{k};
    % Hessian of U = -(x^2+y^2)/2 - (1-mu)/r1 - mu/r2
    d1 = [r(1) + mu; r(2); r(3)];      r1 = norm(d1);
    d2 = [r(1) - 1 + mu; r(2); r(3)];  r2 = norm(d2);
    H = -diag([1 1 0]) + (1-mu)*(eye(3)/r1^3 - 3*(d1*d1')/r1^5) ...
                       +  mu   *(eye(3)/r2^3 - 3*(d2*d2')/r2^5);
    % Jacobian of the first-order system, d(ds)/ds with s = [r; rdot]
    A = [zeros(3), eye(3); -H, Omega];
    lam = eig(A);
    fprintf('\n%s at (%.6f, %+.6f, 0)\nJacobian A =\n', names{k}, r(1), r(2));
    disp(A);
    fprintf('Eigenvalues:\n');
    fprintf('  %+.6f %+.6fi\n', [real(lam) imag(lam)]');
    if max(real(lam)) > 1e-9
        fprintf('Unstable: a deviation grows like exp(%.4f t), i.e., by a factor of e every %.2f days\n', ...
                max(real(lam)), 375190/86400/max(real(lam)));
    else
        fprintf('All eigenvalues are imaginary (max real part %.1e): linearly stable\n', max(real(lam)));
    end
end
