function hg = plotRigidBodyTrajectory(ax, h, J, varargin)
%PLOTRIGIDBODYTRAJECTORY Overlay a rigid body trajectory in momentum space.
%
%   hg = plotRigidBodyTrajectory(ax, h, J) plots the angular momentum
%   trajectory h on the axes ax (typically created by plotMomentumSphere),
%   colored by its energy T = (1/2) h' * (J \ h). Returns a struct of
%   graphics handles.
%
%   ax : target axes ([] for current axes)
%   h  : trajectory, N x 3 (as returned by ode45) or 3 x N
%   J  : 3x3 symmetric positive-definite inertia matrix
%
%   Name-value options:
%     'Color'      line color; 'energy' colors by energy  (default: 'energy')
%     'LineWidth'  line width                              (default: 1.8)
%     'MarkStart'  draw a dot at the initial condition     (default: false)
%
%   Example:
%     J = diag([1 2 3]);
%     ax = plotMomentumSphere(J, 1);
%     rhs  = @(t,h) cross(h, J\h);
%     opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
%     for a = linspace(0.1, pi/2-0.1, 8)
%         [~,h] = ode45(rhs, [0 60], [cos(a); 0; sin(a)], opts);
%         plotRigidBodyTrajectory(ax, h, J);
%     end
%
%   See also PLOTMOMENTUMSPHERE, PLOTENERGYELLIPSOID.

% ---------- Parse inputs ----------
p = inputParser;
p.addParameter('Color', 'energy');
p.addParameter('LineWidth', 1.8, @isnumeric);
p.addParameter('MarkStart', false, @islogical);
p.parse(varargin{:});
opt = p.Results;

if isempty(ax), ax = gca; end
if ~isequal(size(J), [3 3])
    error('plotRigidBodyTrajectory:J', 'J must be a 3x3 inertia matrix.');
end

% Accept N x 3 or 3 x N; treat 3 x 3 as N x 3 (ode45 convention)
if size(h,2) == 3
    h = h.';
elseif size(h,1) ~= 3
    error('plotRigidBodyTrajectory:size', 'h must be N x 3 or 3 x N.');
end

% ---------- Color ----------
h0 = h(:,1);
T0 = 0.5 * (h0.' * (J \ h0));
if ischar(opt.Color) && strcmpi(opt.Color,'energy')
    if isfield(ax.UserData,'rigidBody')
        rb = ax.UserData.rigidBody;  cmap = rb.cmap;  Tlim = rb.Tlim;
    else
        cmap = colormap(ax);
        Jp   = eig((J + J.')/2);
        Tlim = norm(h0)^2 ./ (2*[max(Jp) min(Jp)]);
    end
    nc  = size(cmap,1);
    idx = round(1 + (nc-1)*(T0-Tlim(1))/(Tlim(2)-Tlim(1)));
    c   = cmap(max(1,min(nc,idx)), :);
else
    c = opt.Color;
end

% ---------- Plot ----------
holdState = ishold(ax);
hold(ax,'on');

hg.line = plot3(ax, h(1,:),h(2,:),h(3,:), 'Color',c, ...
    'LineWidth',opt.LineWidth, 'DisplayName',sprintf('T = %.3g', T0));

if opt.MarkStart
    hg.start = plot3(ax, h0(1),h0(2),h0(3), 'o', 'MarkerSize',6, ...
        'MarkerFaceColor',c, 'MarkerEdgeColor','k', 'HandleVisibility','off');
end

axis(ax,'equal');                % keep equal data units on all axes

if ~holdState, hold(ax,'off'); end
end
