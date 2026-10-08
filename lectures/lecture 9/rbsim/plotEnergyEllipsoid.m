function hs = plotEnergyEllipsoid(ax, J, T, varargin)
%PLOTENERGYELLIPSOID Draw rigid body energy ellipsoids in momentum space.
%
%   hs = plotEnergyEllipsoid(ax, J, T) draws the ellipsoid
%       (1/2) h' * (J \ h) = T
%   for each energy level in T on the axes ax, and returns the surface
%   handles. The ellipsoid's axes are the principal axes of J, with
%   semi-axes sqrt(2 T Jp_i), where Jp_i are the principal moments.
%
%   ax : target axes ([] for current axes)
%   J  : 3x3 symmetric positive-definite inertia matrix
%   T  : scalar or vector of energy levels
%
%   Name-value options:
%     'Color'       face color; 'energy' colors by energy  (default: 'energy')
%     'FaceAlpha'   face transparency                     (default: 0.10)
%     'EdgeAlpha'   mesh-line transparency                (default: 0.12)
%     'Resolution'  number of mesh faces around           (default: 60)
%
%   With 'energy' coloring, the colormap and energy limits are taken from
%   the axes if they were set up by plotMomentumSphere.
%
%   Useful energy levels for |h| = hmag, with principal moments
%   Jp1 < Jp2 < Jp3:
%     hmag^2/(2*Jp3)  min energy  (ellipsoid touches sphere at +/- v_3)
%     hmag^2/(2*Jp2)  separatrix
%     hmag^2/(2*Jp1)  max energy  (ellipsoid touches sphere at +/- v_1)
%
%   Example:
%     J = diag([1 2 3]); hmag = 1;
%     ax = plotMomentumSphere(J, hmag);
%     Jp = ax.UserData.rigidBody.Jp;
%     plotEnergyEllipsoid(ax, J, hmag^2/(2*Jp(2)));   % separatrix level
%
%   See also PLOTMOMENTUMSPHERE, PLOTRIGIDBODYTRAJECTORY.

% ---------- Parse inputs ----------
p = inputParser;
p.addParameter('Color', 'energy');
p.addParameter('FaceAlpha', 0.10, @isnumeric);
p.addParameter('EdgeAlpha', 0.12, @isnumeric);
p.addParameter('Resolution', 60, @(x) isnumeric(x) && isscalar(x) && x>=4);
p.parse(varargin{:});
opt = p.Results;

if isempty(ax), ax = gca; end
if ~isequal(size(J), [3 3])
    error('plotEnergyEllipsoid:J', 'J must be a 3x3 inertia matrix.');
end
if any(T(:) <= 0)
    error('plotEnergyEllipsoid:T', 'Energy levels must be positive.');
end

% Principal axes, ascending moments, right-handed
[V, D] = eig((J + J.')/2);
[Jp, ix] = sort(diag(D));
V = V(:, ix);
if det(V) < 0, V(:,3) = -V(:,3); end

% ---------- Color setup ----------
useEnergyColor = ischar(opt.Color) && strcmpi(opt.Color,'energy');
if useEnergyColor
    if isfield(ax.UserData,'rigidBody')
        cmap = ax.UserData.rigidBody.cmap;
        Tlim = ax.UserData.rigidBody.Tlim;
    else
        cmap = colormap(ax);
        Tlim = [min(T(:)) max(T(:))];
        if diff(Tlim) == 0, Tlim = Tlim .* [0.5 1.5]; end
    end
    nc = size(cmap,1);
end

% ---------- Draw ----------
holdState = ishold(ax);
hold(ax,'on');

[xu,yu,zu] = sphere(opt.Resolution);
U  = [xu(:) yu(:) zu(:)].';
sz = size(xu);

T  = T(:).';
hs = gobjects(1, numel(T));
for k = 1:numel(T)
    if useEnergyColor
        idx = round(1 + (nc-1)*(T(k)-Tlim(1))/(Tlim(2)-Tlim(1)));
        c   = cmap(max(1,min(nc,idx)), :);
    else
        c = opt.Color;
    end
    P = V*diag(sqrt(2*T(k)*Jp))*U;
    hs(k) = surf(ax, reshape(P(1,:),sz), reshape(P(2,:),sz), reshape(P(3,:),sz), ...
        'FaceColor',c, 'FaceAlpha',opt.FaceAlpha, ...
        'EdgeColor',c, 'EdgeAlpha',opt.EdgeAlpha, 'FaceLighting','gouraud', ...
        'DisplayName',sprintf('Energy ellipsoid T = %.3g', T(k)));
end

% Make sure the largest ellipsoid fits in the view
r = max(sqrt(2*max(T)*Jp));
lim = axis(ax);
if max(abs(lim)) < 1.1*r
    axis(ax, 1.25*r*[-1 1 -1 1 -1 1]);
end
axis(ax,'equal');                % keep equal data units on all axes

if ~holdState, hold(ax,'off'); end
end
