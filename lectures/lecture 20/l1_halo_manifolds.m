%% Earth-Moon L1 halo orbit and its stable and unstable manifolds (CR3BP)
%
% Dynamics are the same as in 16.07 Problem Set 7, Question 2: rotating
% frame, nondimensional units, and the effective potential
%     U(x,y,z) = -(x^2 + y^2)/2 - (1-mu)/r1 - mu/r2,
% with equations of motion
%     xdd - 2 yd = -dU/dx,   ydd + 2 xd = -dU/dy,   zdd = -dU/dz.
%
% Steps:
%   1. Find the L1 point.
%   2. Compute a periodic halo orbit by differential correction.
%   3. Compute the monodromy matrix and its stable/unstable eigenvectors.
%   4. At N points around the orbit, map the eigenvectors along with the
%      state transition matrix (STM), perturb the state a small distance
%      along them, and propagate: forward in time for the unstable
%      manifold, backward in time for the stable manifold.
%   5. Plot everything in 3D with the Earth and Moon.
%
% All functions are at the end of this file (requires MATLAB R2016b or later).

clear; close all; clc;

%% Parameters
mu    = 0.012150585609624;   % Earth-Moon mass ratio
Lstar = 384400;              % length unit (km)
Tstar = 375190;              % time unit (s)
R_E   = 6378.137 / Lstar;    % Earth radius (nondimensional)
R_M   = 1737.4   / Lstar;    % Moon radius (nondimensional)

z0_halo  = 0.06;             % out-of-plane amplitude at the x-z plane crossing (~23,000 km)
N_points = 40;               % number of points around the halo to seed manifolds from
eps_km   = 50;               % size of the initial perturbation (km)
t_manif  = 2.0 * pi;         % how long to propagate each manifold trajectory (~27 days)

opts_tight = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

%% 1. L1 point: dU/dx = 0 on the x axis between the Earth and the Moon
dUdx_axis = @(x) -x + (1-mu)*(x+mu)./abs(x+mu).^3 + mu*(x-1+mu)./abs(x-1+mu).^3;
x_L1 = fzero(dUdx_axis, [0.5, 1 - mu - 0.01]);
fprintf('L1 is at x = %.10f (%.0f km from the Moon''s center)\n', ...
        x_L1, (1 - mu - x_L1) * Lstar);

%% 2. Halo orbit by differential correction
% The orbit crosses the x-z plane (y = 0) perpendicularly, with xdot = zdot = 0.
% Start at s0 = [x0 0 z0 0 ydot0 0], hold z0 fixed, and adjust x0 and ydot0
% until xdot = zdot = 0 at the next crossing of y = 0 (half a period later).
x0   = 0.8243;    % initial guesses (good for z0 near 0.06; for very
ydot0 = 0.1709;   % different amplitudes, step z0 gradually)

for iter = 1:30
    s0 = [x0; 0; z0_halo; 0; ydot0; 0];
    [t_half, s_half, Phi_half] = propagate_to_y_crossing(s0, mu, opts_tight);
    err = [s_half(4); s_half(6)];                       % want xdot = zdot = 0
    fprintf('  iteration %2d: |[xdot, zdot]| = %.3e\n', iter, norm(err));
    if norm(err) < 1e-11
        break
    end
    f_half = cr3bp_rhs(0, s_half, mu);                  % [xdot ydot zdot xdd ydd zdd]
    % Columns of Phi for the free variables (x0, ydot0), corrected for the
    % change in crossing time needed to stay on y = 0.
    D = [Phi_half(4,1) Phi_half(4,5);
         Phi_half(6,1) Phi_half(6,5)] ...
        - [f_half(4); f_half(6)] * [Phi_half(2,1) Phi_half(2,5)] / f_half(2);
    dv = -D \ err;
    x0    = x0    + dv(1);
    ydot0 = ydot0 + dv(2);
end
T_halo = 2 * t_half;
s0_halo = [x0; 0; z0_halo; 0; ydot0; 0];
fprintf('\nHalo orbit: x0 = %.12f, z0 = %.12f, ydot0 = %.12f\n', x0, z0_halo, ydot0);
fprintf('Period T = %.10f (%.2f days)\n', T_halo, T_halo * Tstar / 86400);

%% 3. Monodromy matrix and stable/unstable directions
% Propagate the state and STM for one full period, saving N_points evenly
% spaced states and STMs along the way.
t_seed = linspace(0, T_halo, N_points + 1);
X0 = [s0_halo; reshape(eye(6), 36, 1)];
[~, X] = ode113(@(t, X) cr3bp_stm_rhs(t, X, mu), t_seed, X0, opts_tight);
Monodromy = reshape(X(end, 7:42), 6, 6);
fprintf('Closure error after one period: %.2e\n', norm(X(end, 1:6)' - s0_halo));

[V, Lam] = eig(Monodromy);
lam = diag(Lam);
[lam_u, i_u] = max(abs(lam));
[lam_s, i_s] = min(abs(lam));
v_u = real(V(:, i_u));
v_s = real(V(:, i_s));
% Fix the signs so that the "+" branch initially moves toward the Moon (+x).
if v_u(1) < 0, v_u = -v_u; end
if v_s(1) < 0, v_s = -v_s; end
fprintf('Monodromy eigenvalues: unstable %.2f, stable %.5f (product %.4f)\n', ...
        lam_u, lam_s, lam_u * lam_s);
fprintf('A small deviation grows by a factor of %.0f every period (%.1f days)\n', ...
        lam_u, T_halo * Tstar / 86400);

%% 4. Manifold trajectories
% At seed point k, the eigenvectors are mapped by the STM from t = 0:
% v_k = Phi(t_k, 0) v. Each is normalized so the position perturbation is
% eps_km, then added (+ branch) or subtracted (- branch) from the state.
eps_nd = eps_km / Lstar;
opts_manif = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
                    'Events', @(t, s) collision_event(t, s, mu, R_E, R_M));

branches = struct('name', {'Unstable, + (Moon side)', 'Unstable, - (Earth side)', ...
                           'Stable, + (Moon side)',   'Stable, - (Earth side)'}, ...
                  'vec',  {'u', 'u', 's', 's'}, ...
                  'sign', {+1, -1, +1, -1}, ...
                  'color', {[0.92 0.41 0.20], [0.93 0.63 0.00], ...   % warm: unstable
                            [0.16 0.47 0.84], [0.11 0.69 0.48]}, ...  % cool: stable
                  'traj', {{}, {}, {}, {}});

for b = 1:numel(branches)
    for k = 1:N_points
        s_k   = X(k, 1:6)';
        Phi_k = reshape(X(k, 7:42), 6, 6);
        if branches(b).vec == 'u'
            v_k = Phi_k * v_u;  t_end =  t_manif;   % unstable: forward in time
        else
            v_k = Phi_k * v_s;  t_end = -t_manif;   % stable: backward in time
        end
        v_k = v_k / norm(v_k(1:3));
        s_pert = s_k + branches(b).sign * eps_nd * v_k;
        [~, s_traj] = ode113(@(t, s) cr3bp_rhs(t, s, mu), [0 t_end], s_pert, opts_manif);
        branches(b).traj{k} = s_traj(:, 1:3);
    end
    fprintf('Propagated %s\n', branches(b).name);
end

%% 5. Plot (in thousands of km, origin at the Earth-Moon barycenter)
sc = Lstar / 1e3;
figure('Color', 'w', 'Position', [100 100 1100 800]); hold on;

h = gobjects(0);
labels = {};
for b = 1:numel(branches)
    for k = 1:N_points
        r = branches(b).traj{k} * sc;
        hl = plot3(r(:,1), r(:,2), r(:,3), '-', 'Color', branches(b).color, 'LineWidth', 0.6);
    end
    h(end+1) = hl;   %#ok<SAGROW>
    labels{end+1} = branches(b).name;   %#ok<SAGROW>
end

% Halo orbit (sample it finely for a smooth curve)
[~, S_halo] = ode113(@(t, s) cr3bp_rhs(t, s, mu), linspace(0, T_halo, 500), s0_halo, opts_tight);
h(end+1) = plot3(S_halo(:,1)*sc, S_halo(:,2)*sc, S_halo(:,3)*sc, 'k-', 'LineWidth', 2.5);
labels{end+1} = 'L1 halo orbit';

% Earth, Moon (true size), and L1
[xs, ys, zs] = sphere(40);
surf((-mu + R_E*xs)*sc, R_E*ys*sc, R_E*zs*sc, ...
     'FaceColor', [0.25 0.45 0.85], 'EdgeColor', 'none');
surf((1 - mu + R_M*xs)*sc, R_M*ys*sc, R_M*zs*sc, ...
     'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
plot3(x_L1*sc, 0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
text(-mu*sc, 0, 25, 'Earth', 'HorizontalAlignment', 'center');
text((1-mu)*sc, 0, 25, 'Moon', 'HorizontalAlignment', 'center');
text(x_L1*sc, 0, -25, 'L_1', 'HorizontalAlignment', 'center');

axis equal; grid on; box on;
xlabel('x (thousand km)'); ylabel('y (thousand km)'); zlabel('z (thousand km)');
title(sprintf(['Earth-Moon L_1 halo orbit (A_z \\approx %.0f km, T = %.1f days) ' ...
               'and its invariant manifolds'], z0_halo*Lstar, T_halo*Tstar/86400));
legend(h, labels, 'Location', 'northeastoutside');
lighting gouraud; camlight;
view(-35, 25);

%% ========================================================================
%  Local functions
%  ========================================================================

function ds = cr3bp_rhs(~, s, mu)
% CR3BP equations of motion in the rotating frame (nondimensional).
    gU = grad_U(s(1:3), mu);
    ds = [s(4); s(5); s(6);
          2*s(5) - gU(1);
         -2*s(4) - gU(2);
                 - gU(3)];
end

function g = grad_U(r, mu)
% Gradient of U = -(x^2 + y^2)/2 - (1-mu)/r1 - mu/r2.
    x = r(1); y = r(2); z = r(3);
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    g = [-x + (1-mu)*(x + mu)/r1^3 + mu*(x - 1 + mu)/r2^3;
         -y + (1-mu)*y/r1^3        + mu*y/r2^3;
              (1-mu)*z/r1^3        + mu*z/r2^3];
end

function H = hess_U(r, mu)
% Hessian of U (3x3, symmetric).
    d1 = [r(1) + mu; r(2); r(3)];      r1 = norm(d1);
    d2 = [r(1) - 1 + mu; r(2); r(3)];  r2 = norm(d2);
    H = -diag([1 1 0]) ...
        + (1-mu) * (eye(3)/r1^3 - 3*(d1*d1')/r1^5) ...
        +  mu    * (eye(3)/r2^3 - 3*(d2*d2')/r2^5);
end

function dX = cr3bp_stm_rhs(t, X, mu)
% State plus state transition matrix: dPhi/dt = A(t) Phi.
    s   = X(1:6);
    Phi = reshape(X(7:42), 6, 6);
    A = [zeros(3), eye(3);
         -hess_U(s(1:3), mu), [0 2 0; -2 0 0; 0 0 0]];
    dX = [cr3bp_rhs(t, s, mu); reshape(A * Phi, 36, 1)];
end

function [t_c, s_c, Phi_c] = propagate_to_y_crossing(s0, mu, opts)
% Propagate state + STM from s0 (on y = 0, ydot > 0) to the next crossing
% of y = 0 going downward, which is half a period for a symmetric orbit.
    opts = odeset(opts, 'Events', @y_crossing_event);
    X0 = [s0; reshape(eye(6), 36, 1)];
    [~, ~, t_e, X_e] = ode113(@(t, X) cr3bp_stm_rhs(t, X, mu), [0 10], X0, opts);
    t_c   = t_e(end);
    s_c   = X_e(end, 1:6)';
    Phi_c = reshape(X_e(end, 7:42), 6, 6);
end

function [value, isterminal, direction] = y_crossing_event(~, X)
    value = X(2);         % y = 0
    isterminal = 1;
    direction = -1;       % y decreasing through zero
end

function [value, isterminal, direction] = collision_event(~, s, mu, R_E, R_M)
% Stop a manifold trajectory if it hits the Earth or the Moon.
    r1 = norm([s(1) + mu; s(2); s(3)]);
    r2 = norm([s(1) - 1 + mu; s(2); s(3)]);
    value = [r1 - R_E; r2 - R_M];
    isterminal = [1; 1];
    direction = [-1; -1];
end
