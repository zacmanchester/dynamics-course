%% Effective potential and zero-velocity (Hill) curves of the Earth-Moon CR3BP
%
% Same dynamics and sign convention as Problem Set 7: in the rotating frame,
%     U(x,y,z) = -(x^2 + y^2)/2 - (1-mu)/r1 - mu/r2,
% and the energy E = (xd^2 + yd^2 + zd^2)/2 + U is conserved.
% Since the kinetic energy can't be negative, a spacecraft with energy E can
% only be where U(x,y,z) <= E. The surface U(x,y,z) = E, where the speed in
% the rotating frame would be zero, is the zero-velocity surface (Hill
% surface). It bounds the "forbidden" region U > E.
%
% Figure 1 plots U over the x-y plane as a surface, with the energy level E
% as a plane. Figure 2 shows the zero-velocity curve in the x-y plane, with
% the forbidden region shaded.
%
% Try changing E_frac below to see the allowed regions change:
%   E < E_L1        : the spacecraft is trapped near the Earth or near the Moon
%   E_L1 < E < E_L2 : a "neck" opens at L1, connecting the Earth and Moon regions
%   E_L2 < E < E_L3 : a second neck opens at L2, so it can escape past the Moon
%   E_L3 < E < E_L4 : the necks merge; only regions around L4 and L5 are forbidden
%   E > E_L4        : no forbidden region in the x-y plane

clear; close all;

mu = 0.012150585609624;          % Earth-Moon mass ratio
L_star = 384400;                 % length unit (km), for labeling

U = @(x, y, z) -0.5*(x.^2 + y.^2) ...
               - (1-mu)./sqrt((x+mu).^2 + y.^2 + z.^2) ...
               - mu./sqrt((x-1+mu).^2 + y.^2 + z.^2);

%% Lagrange points
% Collinear points: solve dU/dx = 0 on the x axis (see lagrange_point_eigs.m).
% Each bracket must contain exactly one sign change of dU/dx.
dUdx_axis = @(x) -x + (1-mu)*(x+mu)./abs(x+mu).^3 + mu*(x-1+mu)./abs(x-1+mu).^3;
xL = zeros(5, 1);  yL = zeros(5, 1);
xL(1) = fzero(dUdx_axis, [0.5, 1 - mu - 0.01]);     % between the Earth and Moon
xL(2) = fzero(dUdx_axis, [1 - mu + 0.01, 1.5]);     % beyond the Moon
xL(3) = fzero(dUdx_axis, [-1.5, -mu - 0.01]);       % opposite side of the Earth
xL(4) = 0.5 - mu;  yL(4) =  sqrt(3)/2;              % equilateral triangle points
xL(5) = 0.5 - mu;  yL(5) = -sqrt(3)/2;
E_L = U(xL, yL, 0*xL);                              % energy of each point
for k = 1:5
    fprintf('L%d at (%+.6f, %+.6f):  E = U(L%d) = %.6f\n', k, xL(k), yL(k), k, E_L(k));
end

%% Choose the energy level
% E_frac = 0 puts E at E_L1, and E_frac = 1 puts it at E_L2.
E_frac = 1.5; %-5.0; %0.5; 
E = E_L(1) + E_frac * (E_L(2) - E_L(1));
fprintf('\nEnergy level: E = %.6f\n', E);

%% Effective potential on a grid in the x-y plane (z = 0)
n2 = 601;
[X2, Y2] = meshgrid(linspace(-1.5, 1.5, n2) + 1e-4, linspace(-1.5, 1.5, n2) + 1e-4);
U2 = U(X2, Y2, 0*X2);

%% Effective potential over the x-y plane (z = 0), plotted as a surface
% The height of the surface is U(x,y,0). U goes to -infinity at the Earth
% and Moon, so we clip it at U_floor to keep the plot readable. The
% translucent plane is the energy level E: the spacecraft can only be where
% the surface lies below the plane, and the black curve where they meet is
% the zero-velocity curve from the 2D plot. L4 and L5 are hilltops of U,
% and L1, L2, and L3 are saddle points.
U_floor = -2.0;
U2c = max(U2, U_floor);

figure('Color', 'w', 'Position', [100 100 900 700]); hold on;
surf(X2, Y2, U2c, 'EdgeColor', 'none', 'FaceAlpha', 0.95);
% Energy plane and zero-velocity curve
surf([-1.5 1.5; -1.5 1.5], [-1.5 -1.5; 1.5 1.5], E*ones(2), ...
     'FaceColor', [0.85 0.33 0.10], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
contour3(X2, Y2, U2, [E E], 'k', 'LineWidth', 1.5);
% Lagrange points at their potential values
plot3(xL, yL, E_L, 'k.', 'MarkerSize', 18);
for k = 1:5
    text(xL(k), yL(k), E_L(k) + 0.04, sprintf('L_%d', k), ...
         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
text(-mu, 0, U_floor + 0.03, 'Earth', 'HorizontalAlignment', 'center');
text(1 - mu, 0, U_floor + 0.03, 'Moon', 'HorizontalAlignment', 'center');

xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([U_floor, -1.4]);
xlabel('x'); ylabel('y'); zlabel('U(x, y, 0)');
title(sprintf('Effective potential in the x-y plane (clipped at %.1f); plane: E = %.4f', U_floor, E));
cb = colorbar; ylabel(cb, 'U');
grid on; box on;
view(-35, 40);

%% Slice in the x-y plane (z = 0)
% Shade the forbidden region U > E and draw the zero-velocity curve.
figure('Color', 'w', 'Position', [150 150 700 650]); hold on;
contourf(X2, Y2, double(U2 > E), [0.5 0.5], 'LineStyle', 'none');
colormap([1 1 1; 0.75 0.82 0.92]);           % allowed: white, forbidden: light blue
contour(X2, Y2, U2, [E E], 'LineColor', [0.16 0.47 0.84], 'LineWidth', 1.5);
plot(-mu, 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.25 0.45 0.85], 'MarkerEdgeColor', 'k');
plot(1 - mu, 0, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k');
plot(xL, yL, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
for k = 1:5
    text(xL(k) + 0.04, yL(k) + 0.06, sprintf('L_%d', k));
end
text(-mu, -0.1, 'Earth', 'HorizontalAlignment', 'center');
text(1 - mu, -0.1, 'Moon', 'HorizontalAlignment', 'center');
axis equal; box on;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
xlabel('x'); ylabel('y');
title(sprintf('Zero-velocity curve in the x-y plane, E = %.4f (shaded: forbidden)', E));
