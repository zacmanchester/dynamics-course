%% A chaotic trajectory in the Earth-Moon CR3BP
% Same dynamics and sign convention as Problem Set 7. The initial condition
% starts near the Moon with an energy between E_L1 and E_L2: the neck at L1
% is open but the one at L2 is closed. The spacecraft can't escape, so it
% wanders between the Earth and Moon regions, passing through the L1 neck at
% irregular times. A second trajectory starting 1 km away shows the
% sensitive dependence on initial conditions.

mu = 0.012150585609624;  L_star = 384400;  T_star_days = 375190/86400;
gradU = @(x, y) [-x + (1-mu)*(x+mu)/((x+mu)^2+y^2)^1.5 + mu*(x-1+mu)/((x-1+mu)^2+y^2)^1.5;
                 -y + (1-mu)*y/((x+mu)^2+y^2)^1.5      + mu*y/((x-1+mu)^2+y^2)^1.5];
rhs = @(t, s) [s(3); s(4); [2*s(4); -2*s(3)] - gradU(s(1), s(2))];   % planar EOM
U = @(x, y) -0.5*(x.^2 + y.^2) - (1-mu)./sqrt((x+mu).^2 + y.^2) - mu./sqrt((x-1+mu).^2 + y.^2);

% Initial condition [x0 y0 xdot0 ydot0] (planar, z = zdot = 0)
s0 = [0.90; 0; 0; -0.2719726169];
fprintf('Energy E = %.6f  (E_L1 = -1.594171, E_L2 = -1.586080)\n', 0.5*s0(4)^2 + U(s0(1), s0(2)));

tf = 100;                                    % about 434 days
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
t = linspace(0, tf, 20000);
[~, S1] = ode45(rhs, t, s0, opts);
[~, S2] = ode45(rhs, t, s0 + [1/L_star; 0; 0; 0], opts);   % start 1 km away in x

% Trajectories in the rotating frame
figure('Color', 'w', 'Position', [100 100 1100 480]);
subplot(1, 2, 1); hold on;
plot(S1(:,1), S1(:,2), 'Color', [0.16 0.47 0.84], 'LineWidth', 0.7);
plot(S2(:,1), S2(:,2), 'Color', [0.92 0.41 0.20], 'LineWidth', 0.7);
[Xg, Yg] = meshgrid(linspace(-1.3, 1.3, 400));
contour(Xg, Yg, U(Xg, Yg), [1 1] * (0.5*s0(4)^2 + U(s0(1), s0(2))), 'k');  % zero-velocity curve
plot(-mu, 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.25 0.45 0.85], 'MarkerEdgeColor', 'k');
plot(1-mu, 0, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k');
plot(0.836915, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);  text(0.85, 0.06, 'L_1');
axis equal; box on; xlim([-1.3 1.3]); ylim([-1.3 1.3]);
xlabel('x'); ylabel('y');
legend('Nominal', 'Started 1 km away', 'Zero-velocity curve', 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Trajectories in the rotating frame');

% Separation between the two trajectories
subplot(1, 2, 2);
d = sqrt((S1(:,1) - S2(:,1)).^2 + (S1(:,2) - S2(:,2)).^2) * L_star;
semilogy(t * T_star_days, d, 'k', 'LineWidth', 1);
grid on; xlabel('Time (days)'); ylabel('Separation (km)');
title('Distance between the two trajectories');
