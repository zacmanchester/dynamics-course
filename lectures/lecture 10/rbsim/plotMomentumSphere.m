function ax = plotMomentumSphere(J, hmag, varargin)
%PLOTMOMENTUMSPHERE Draw the momentum sphere of a free rigid body.
%
%   ax = plotMomentumSphere(J, hmag) draws, in body angular momentum space,
%   the momentum sphere |h| = hmag together with the separatrices and the
%   six equilibria (the principal axes), and returns the axes handle.
%   Energy ellipsoids are drawn separately with plotEnergyEllipsoid.
%
%   J    : 3x3 symmetric positive-definite inertia matrix with distinct
%          principal moments. J need not be diagonal: the equilibria and
%          separatrices are aligned with its eigenvectors.
%   hmag : magnitude of the angular momentum |h|
%
%   Dynamics (body frame):  dh/dt = h x (J \ h)
%   Conserved:  |h|^2  and  T = (1/2) h' * (J \ h)
%
%   Name-value options:
%     'Parent'          axes to draw into              (default: new figure)
%     'ShowSeparatrix'  true/false                     (default: true)
%     'ShowEquilibria'  true/false                     (default: true)
%     'SphereAlpha'     sphere face transparency       (default: 0.25)
%     'Colormap'        colormap used for energy       (default: parula(256))
%
%   Equilibria are labeled v_1, v_2, v_3: the principal axes ordered by
%   increasing principal moment (v_2 is the unstable intermediate axis).
%
%   The energy range [hmag^2/(2 Jmax), hmag^2/(2 Jmin)] is set as the axes
%   color limits and stored in ax.UserData.rigidBody, so ellipsoids and
%   trajectories added later are colored consistently.
%
%   Example:
%     J = diag([1 2 3]); hmag = 1;
%     ax = plotMomentumSphere(J, hmag);
%     plotEnergyEllipsoid(ax, J, [0.2 0.25 0.35]);
%     [~,h] = ode45(@(t,h) cross(h, J\h), [0 50], [0.9; 0; sqrt(1-0.81)]);
%     plotRigidBodyTrajectory(ax, h, J);
%
%   See also PLOTENERGYELLIPSOID, PLOTRIGIDBODYTRAJECTORY.

% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('J', @(x) isnumeric(x) && isequal(size(x),[3 3]));
p.addRequired('hmag', @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('Parent', [], @(x) isempty(x) || isgraphics(x,'axes'));
p.addParameter('ShowSeparatrix', true, @islogical);
p.addParameter('ShowEquilibria', true, @islogical);
p.addParameter('SphereAlpha', 0.25, @isnumeric);
p.addParameter('Colormap', parula(256), @(x) isnumeric(x) && size(x,2)==3);
p.parse(J, hmag, varargin{:});
opt = p.Results;

[V, Jp] = principalAxes(J);      % J = V*diag(Jp)*V', Jp ascending

Tlo = hmag^2/(2*Jp(3));          % lowest energy  (max-inertia axis)
Thi = hmag^2/(2*Jp(1));          % highest energy (min-inertia axis)

% ---------- Axes ----------
if isempty(opt.Parent)
    figure('Color','w');
    ax = axes;
else
    ax = opt.Parent;
end
holdState = ishold(ax);
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
xlabel(ax,'h_1'); ylabel(ax,'h_2'); zlabel(ax,'h_3');

% ---------- Momentum sphere ----------
[xs,ys,zs] = sphere(80);
surf(ax, hmag*xs, hmag*ys, hmag*zs, 'FaceColor',[0.82 0.84 0.92], ...
    'FaceAlpha',opt.SphereAlpha, 'EdgeColor','none', ...
    'FaceLighting','gouraud', 'DisplayName','Momentum sphere');

% ---------- Separatrices ----------
% Two great circles through +/- v_2, in the planes spanned by v_2 and
%   cos(aSep) v_1 +/- sin(aSep) v_3,
%   tan^2(aSep) = (1/Jp1 - 1/Jp2) / (1/Jp2 - 1/Jp3)
if opt.ShowSeparatrix
    aSep = atan(sqrt((1/Jp(1)-1/Jp(2)) / (1/Jp(2)-1/Jp(3))));
    s = linspace(0, 2*pi, 400);
    for sgn = [-1 1]
        u = cos(aSep)*V(:,1) + sgn*sin(aSep)*V(:,3);
        C = hmag*(u*cos(s) + V(:,2)*sin(s));
        hl = plot3(ax, C(1,:),C(2,:),C(3,:), 'k--', 'LineWidth',1.6, ...
            'DisplayName','Separatrix');
        if sgn == 1, hl.HandleVisibility = 'off'; end
    end
end

% ---------- Equilibria (principal axes) ----------
if opt.ShowEquilibria
    E = hmag*[V, -V];                    % +v1 +v2 +v3 -v1 -v2 -v3
    stable = logical([1 0 1 1 0 1]);
    plot3(ax, E(1,stable),E(2,stable),E(3,stable), 'o', 'MarkerSize',10, ...
        'MarkerFaceColor',[0.1 0.7 0.2], 'MarkerEdgeColor','k', ...
        'DisplayName','Stable equilibria');
    plot3(ax, E(1,~stable),E(2,~stable),E(3,~stable), 'o', 'MarkerSize',10, ...
        'MarkerFaceColor',[0.85 0.1 0.1], 'MarkerEdgeColor','k', ...
        'DisplayName','Unstable equilibria');
    labels = {'+v_1','+v_2','+v_3','-v_1','-v_2','-v_3'};
    for k = 1:6
        text(ax, 1.12*E(1,k),1.12*E(2,k),1.12*E(3,k), labels{k}, ...
            'FontSize',11, 'FontWeight','bold', 'HorizontalAlignment','center');
    end
end

% ---------- Finishing ----------
colormap(ax, opt.Colormap);
caxis(ax, [Tlo Thi]);
cb = colorbar(ax); cb.Label.String = 'Energy T';
lim = 1.25*hmag;
axis(ax, lim*[-1 1 -1 1 -1 1]);
axis(ax,'equal');                % equal data units on all axes (after limits)
view(ax, [135 22]);
lighting(ax,'gouraud');
if isempty(findobj(ax,'Type','light')), camlight(ax,'headlight'); end
title(ax, sprintf('Principal moments [%g, %g, %g],  |h| = %g', Jp, hmag));

% Store parameters for plotEnergyEllipsoid / plotRigidBodyTrajectory
ax.UserData.rigidBody = struct('J',J, 'V',V, 'Jp',Jp, 'hmag',hmag, ...
                               'Tlim',[Tlo Thi], 'cmap',opt.Colormap);

if ~holdState, hold(ax,'off'); end
end

% =====================================================================
function [V, Jp] = principalAxes(J)
%PRINCIPALAXES Validate J and return right-handed principal axes.
if norm(J - J.', 'fro') > 1e-10*norm(J, 'fro')
    error('plotMomentumSphere:notSymmetric', 'J must be symmetric.');
end
[V, D] = eig((J + J.')/2);
[Jp, ix] = sort(diag(D));
V = V(:, ix);
if Jp(1) <= 0
    error('plotMomentumSphere:notPosDef', 'J must be positive definite.');
end
if min(diff(Jp)) < 1e-8*Jp(3)
    error('plotMomentumSphere:degenerate', ...
          'Principal moments must be distinct (asymmetric body).');
end
if det(V) < 0, V(:,3) = -V(:,3); end     % make the frame right-handed
end
