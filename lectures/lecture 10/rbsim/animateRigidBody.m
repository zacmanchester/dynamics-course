function [ax, axS] = animateRigidBody(t, R, body, varargin)
%ANIMATERIGIDBODY Animate a tumbling rigid body in the inertial frame.
%
%   ax = animateRigidBody(t, R, body) animates the body described by the
%   struct body (e.g. from makeTBody) as it rotates through the attitudes
%   R(:,:,k) at times t(k). Returns the axes handle.
%
%   [ax, axS] = animateRigidBody(..., 'hBody', hB, 'J', J) also animates
%   the body-frame angular momentum on the momentum sphere in a second
%   panel (axS), in sync with the body. This panel uses plotMomentumSphere
%   and plotRigidBodyTrajectory, which must be on the MATLAB path.
%
%   t    : N-vector of times (uniform spacing recommended)
%   R    : 3x3xN rotation matrices mapping body-frame vectors to the
%          inertial frame (v_inertial = R(:,:,k) * v_body)
%   body : struct with fields
%            Vertices  Mx3 body-frame vertices (COM at origin)
%            Faces     patch faces
%          and optionally
%            FacePart  part index per face (each part gets its own color)
%            tip       3x1 body point to trace with a trail
%
%   Name-value options:
%     'h'           3x1 inertial angular momentum, drawn as a fixed arrow
%     'omega'       3xN inertial angular velocity, drawn as a moving arrow
%     'Speed'       playback speed relative to real time   (default: 1)
%     'ShowAxes'    draw the body-fixed x, y, z axes         (default: true)
%     'Trail'       trace the path of body.tip               (default: true)
%     'TrailLength' number of past frames kept in the trail  (default: Inf)
%     'PartColors'  K x 3 colors for the body parts
%     'hBody'       3xN body-frame angular momentum; adds the momentum
%                   sphere panel (requires 'J')
%     'J'           3x3 inertia matrix (for the momentum sphere panel)
%     'Energy'      1xN energy used to color the sphere trace point by
%                   point (default: (1/2) h' * (J \ h) along hBody).
%                   Useful for dissipative systems, where energy decays.
%     'ShowEnergy'  add an energy-vs-time panel below      (default: false)
%                   (uses 'Energy', or 'hBody' and 'J')
%     'Parent'      axes for the body                        (default: new figure)
%     'SphereAxes'  axes for the momentum sphere (required with 'hBody'
%                   if 'Parent' is given)
%     'EnergyAxes'  axes for the energy plot (required with 'ShowEnergy'
%                   if 'Parent' is given)
%     'VideoFile'   save an MPEG-4 video to this file        (default: '')
%
%   On screen, playback is paced to real time (scaled by 'Speed'), and
%   frames are skipped if rendering falls behind. When 'VideoFile' is set,
%   every frame is written and the video frame rate is Speed/dt.
%
%   Example:
%     body = makeTBody();
%     t = linspace(0, 5, 300);
%     R = zeros(3,3,numel(t));
%     for k = 1:numel(t)
%         a = 2*pi*t(k);
%         R(:,:,k) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
%     end
%     animateRigidBody(t, R, body);
%
%   See also MAKETBODY, PLOTMOMENTUMSPHERE, PLOTRIGIDBODYTRAJECTORY.

% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('t', @(x) isnumeric(x) && isvector(x));
p.addRequired('R', @(x) isnumeric(x) && size(x,1)==3 && size(x,2)==3);
p.addRequired('body', @(x) isstruct(x) && isfield(x,'Vertices') && isfield(x,'Faces'));
p.addParameter('h', [], @(x) isempty(x) || numel(x)==3);
p.addParameter('omega', [], @(x) isempty(x) || size(x,1)==3);
p.addParameter('Speed', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ShowAxes', true, @islogical);
p.addParameter('Trail', true, @islogical);
p.addParameter('TrailLength', Inf, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('PartColors', [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19], ...
               @(x) isnumeric(x) && size(x,2)==3);
p.addParameter('hBody', [], @(x) isempty(x) || size(x,1)==3);
p.addParameter('J', [], @(x) isempty(x) || isequal(size(x),[3 3]));
p.addParameter('Energy', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('Parent', [], @(x) isempty(x) || isgraphics(x,'axes'));
p.addParameter('SphereAxes', [], @(x) isempty(x) || isgraphics(x,'axes'));
p.addParameter('ShowEnergy', false, @islogical);
p.addParameter('EnergyAxes', [], @(x) isempty(x) || isgraphics(x,'axes'));
p.addParameter('VideoFile', '', @(x) ischar(x) || isstring(x));
p.parse(t, R, body, varargin{:});
opt = p.Results;

t = t(:);
N = numel(t);
if size(R,3) ~= N
    error('animateRigidBody:size', 'R must be 3x3xN with N = numel(t).');
end
if ~isempty(opt.omega) && size(opt.omega,2) ~= N
    error('animateRigidBody:size', 'omega must be 3xN with N = numel(t).');
end
showSphere = ~isempty(opt.hBody);
if showSphere
    if size(opt.hBody,2) ~= N
        error('animateRigidBody:size', 'hBody must be 3xN with N = numel(t).');
    end
    if isempty(opt.J)
        error('animateRigidBody:J', 'The momentum sphere panel needs ''J''.');
    end
    if ~isempty(opt.Parent) && isempty(opt.SphereAxes)
        error('animateRigidBody:axes', ...
              'When ''Parent'' is given, also give ''SphereAxes'' for the momentum sphere.');
    end
end
showEnergy = opt.ShowEnergy;
if ~isempty(opt.Energy) && numel(opt.Energy) ~= N
    error('animateRigidBody:size', 'Energy must have N = numel(t) entries.');
end
% Energy at every sample (colors the sphere trace; drawn in the energy panel)
Tk = [];
if ~isempty(opt.Energy)
    Tk = opt.Energy(:).';
elseif ~isempty(opt.hBody) && ~isempty(opt.J)
    Tk = 0.5*sum(opt.hBody .* (opt.J\opt.hBody), 1);
end
if showEnergy
    if isempty(Tk)
        error('animateRigidBody:energy', ...
              'The energy panel needs ''Energy'', or ''hBody'' and ''J''.');
    end
    if ~isempty(opt.Parent) && isempty(opt.EnergyAxes)
        error('animateRigidBody:axes', ...
              'When ''Parent'' is given, also give ''EnergyAxes'' for the energy plot.');
    end
end
hasTip = isfield(body,'tip') && ~isempty(body.tip);
doTrail = opt.Trail && hasTip;

% ---------- Scene scale ----------
Vb   = body.Vertices;                     % M x 3, body frame
Lmax = max(vecnorm(Vb, 2, 2));
lim  = 1.6*Lmax;
arrowLen = 1.3*Lmax;

% ---------- Axes ----------
axS = gobjects(0);  axE = gobjects(0);
if isempty(opt.Parent)
    if showSphere && showEnergy
        fig = figure('Color','w','Position',[80 60 1450 900]);
        ax  = axes(fig, 'Position',[0.04 0.33 0.43 0.60]);
        axS = axes(fig, 'Position',[0.53 0.33 0.40 0.60]);
        axE = axes(fig, 'Position',[0.07 0.07 0.86 0.18]);
    elseif showSphere
        fig = figure('Color','w','Position',[80 120 1450 680]);
        ax  = subplot(1,2,1,'Parent',fig);
        axS = subplot(1,2,2,'Parent',fig);
    elseif showEnergy
        fig = figure('Color','w','Position',[150 60 800 950]);
        ax  = axes(fig, 'Position',[0.10 0.33 0.80 0.60]);
        axE = axes(fig, 'Position',[0.12 0.07 0.80 0.18]);
    else
        fig = figure('Color','w','Position',[150 150 800 750]);
        ax  = axes(fig);
    end
else
    ax  = opt.Parent;
    fig = ancestor(ax,'figure');
    if showSphere, axS = opt.SphereAxes; end
    if showEnergy, axE = opt.EnergyAxes; end
end
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
axis(ax, lim*[-1 1 -1 1 -1 1]);
axis(ax,'equal');
xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
view(ax, [125 20]);

% ---------- Body patch (per-face colors by part) ----------
nF = size(body.Faces,1);
if isfield(body,'FacePart') && numel(body.FacePart) == nF
    part = body.FacePart(:);
else
    part = ones(nF,1);
end
cols = opt.PartColors(mod(part-1, size(opt.PartColors,1)) + 1, :);

hBody = patch(ax, 'Vertices', (R(:,:,1)*Vb.').', 'Faces', body.Faces, ...
    'FaceVertexCData', cols, 'FaceColor','flat', 'EdgeColor','k', ...
    'LineWidth',0.75, 'FaceLighting','gouraud', 'AmbientStrength',0.4);
camlight(ax,'headlight'); lighting(ax,'gouraud');

% ---------- Body axes ----------
axCols = [0.8 0 0; 0 0.6 0; 0 0 0.8];
hAx = gobjects(1,3);
if opt.ShowAxes
    for i = 1:3
        e = R(:,i,1)*0.9*Lmax;
        hAx(i) = plot3(ax, [0 e(1)], [0 e(2)], [0 e(3)], '-', ...
            'Color',axCols(i,:), 'LineWidth',2);
    end
end

% ---------- Angular momentum (fixed) and angular velocity (moving) ----------
if ~isempty(opt.h)
    hv = opt.h(:)/norm(opt.h)*arrowLen;
    quiver3(ax, 0,0,0, hv(1),hv(2),hv(3), 0, 'Color','k', ...
        'LineWidth',2.5, 'MaxHeadSize',0.4);
    text(ax, 1.08*hv(1), 1.08*hv(2), 1.08*hv(3), 'h', ...
        'FontSize',13, 'FontWeight','bold');
end
hW = gobjects(1);
if ~isempty(opt.omega)
    wScale = arrowLen/max(vecnorm(opt.omega, 2, 1));
    w1 = opt.omega(:,1)*wScale;
    hW = quiver3(ax, 0,0,0, w1(1),w1(2),w1(3), 0, 'Color',[0.6 0 0.6], ...
        'LineWidth',2, 'LineStyle','-', 'MaxHeadSize',0.4);
end

% ---------- Trail of the tip ----------
if doTrail
    tipPath = zeros(3,N);
    for k = 1:N, tipPath(:,k) = R(:,:,k)*body.tip(:); end
    hTrail = plot3(ax, tipPath(1,1), tipPath(2,1), tipPath(3,1), '-', ...
        'Color',[0.3 0.3 0.3 0.6], 'LineWidth',1.2);
end

hTitle = title(ax, sprintf('t = %.2f', t(1)));

% ---------- Momentum sphere panel (body frame) ----------
if showSphere
    hB   = opt.hBody;
    hmag = norm(hB(:,1));
    plotMomentumSphere(opt.J, hmag, 'Parent', axS);
    % Full path, faint, for reference
    plotRigidBodyTrajectory(axS, hB, opt.J, 'Color',[0.55 0.55 0.55], 'LineWidth',0.8);
    rb   = axS.UserData.rigidBody;
    nc   = size(rb.cmap,1);
    Tcol = @(T) rb.cmap(max(1,min(nc, round(1 + (nc-1)*(T-rb.Tlim(1))/diff(rb.Tlim)))), :);
    hold(axS,'on');
    % Trace drawn as a degenerate surface so its color can vary along it
    hSTrace = surface(axS, [hB(1,1);hB(1,1)], [hB(2,1);hB(2,1)], [hB(3,1);hB(3,1)], ...
                      [Tk(1);Tk(1)], 'FaceColor','none', 'EdgeColor','interp', ...
                      'LineWidth',2.5);
    hSVec   = plot3(axS, [0 hB(1,1)],[0 hB(2,1)],[0 hB(3,1)], '-', ...
                    'Color','k', 'LineWidth',1.5);
    hSDot   = plot3(axS, hB(1,1),hB(2,1),hB(3,1), 'o', 'MarkerSize',9, ...
                    'MarkerFaceColor',Tcol(Tk(1)), 'MarkerEdgeColor','k', 'LineWidth',1.2);
    title(axS, 'Body-frame angular momentum h');
    xlabel(axS,'h_x'); ylabel(axS,'h_y'); zlabel(axS,'h_z');
end

% ---------- Energy vs time panel ----------
if showEnergy
    hold(axE,'on'); grid(axE,'on'); box(axE,'on');
    eCol = [0.00 0.45 0.74];
    plot(axE, t, Tk, '-', 'Color',[0.75 0.75 0.75], 'LineWidth',1);   % full history, faint
    Tlo = []; Tsep = []; Thi = [];
    if ~isempty(opt.hBody) && ~isempty(opt.J)
        % Reference levels from the (locked) inertia J at this |h|
        Jp   = sort(eig((opt.J + opt.J.')/2));
        hm2  = norm(opt.hBody(:,1))^2;
        Tlo  = hm2/(2*Jp(3));  Tsep = hm2/(2*Jp(2));  Thi = hm2/(2*Jp(1));
    end
    yl = [min([Tk, Tlo]), max(Tk)];
    if diff(yl) < 1e-9*max(abs(yl)), yl = yl(1)*[0.9 1.1] + [-eps eps]; end
    yl = yl + 0.08*diff(yl)*[-1 1];
    refT = [Tlo, Tsep, Thi];
    refN = {'T_{min}', 'T_{sep}', 'T_{max}'};
    if isempty(Tlo), refN = {}; end
    for i = 1:numel(refT)
        if refT(i) >= yl(1) && refT(i) <= yl(2)
            plot(axE, t([1 end]), refT(i)*[1 1], 'k--', 'LineWidth',0.8);
            text(axE, t(end), refT(i), ['  ' refN{i}], 'FontSize',9, ...
                 'VerticalAlignment','middle', 'Clipping','off');
        end
    end
    hETrace = plot(axE, t(1), Tk(1), '-', 'Color',eCol, 'LineWidth',2);
    hEDot   = plot(axE, t(1), Tk(1), 'o', 'MarkerSize',7, ...
                   'MarkerFaceColor',eCol, 'MarkerEdgeColor','k');
    xlim(axE, [t(1) t(end)]);  ylim(axE, yl);
    xlabel(axE, 't');  ylabel(axE, 'Energy T');
    title(axE, 'Kinetic energy');
end

% ---------- Legend-free annotation ----------
info = {};
if opt.ShowAxes, info{end+1} = 'red/green/blue: body x/y/z axes'; end
if ~isempty(opt.h), info{end+1} = 'black: angular momentum h (fixed)'; end
if ~isempty(opt.omega), info{end+1} = 'purple: angular velocity \omega'; end
if showSphere, info{end+1} = 'right: body-frame h on the momentum sphere'; end
if ~isempty(info)
    if showEnergy, infoPos = [0.01 0.93 0.5 0.07]; else, infoPos = [0.01 0.01 0.6 0.08]; end
    annotation(fig, 'textbox', infoPos, 'String', info, ...
        'EdgeColor','none', 'FontSize',9);
end

% ---------- Frame update ----------
    function drawFrame(k)
        Rk = R(:,:,k);
        hBody.Vertices = (Rk*Vb.').';
        if opt.ShowAxes
            for ii = 1:3
                e = Rk(:,ii)*0.9*Lmax;
                set(hAx(ii), 'XData',[0 e(1)], 'YData',[0 e(2)], 'ZData',[0 e(3)]);
            end
        end
        if ~isempty(opt.omega)
            wk = opt.omega(:,k)*wScale;
            set(hW, 'UData',wk(1), 'VData',wk(2), 'WData',wk(3));
        end
        if doTrail
            k0 = max(1, k - opt.TrailLength + 1);
            set(hTrail, 'XData',tipPath(1,k0:k), 'YData',tipPath(2,k0:k), ...
                        'ZData',tipPath(3,k0:k));
        end
        hTitle.String = sprintf('t = %.2f', t(k));
        if showEnergy
            set(hETrace, 'XData',t(1:k), 'YData',Tk(1:k));
            set(hEDot,   'XData',t(k),   'YData',Tk(k));
        end
        if showSphere
            set(hSTrace, 'XData',[hB(1,1:k);hB(1,1:k)], 'YData',[hB(2,1:k);hB(2,1:k)], ...
                         'ZData',[hB(3,1:k);hB(3,1:k)], 'CData',[Tk(1:k);Tk(1:k)]);
            set(hSVec,   'XData',[0 hB(1,k)], 'YData',[0 hB(2,k)], 'ZData',[0 hB(3,k)]);
            set(hSDot,   'XData',hB(1,k), 'YData',hB(2,k), 'ZData',hB(3,k), ...
                         'MarkerFaceColor',Tcol(Tk(k)));
        end
    end

% ---------- Lock equal axes for the whole animation ----------
% Re-apply after all graphics are created so nothing added above (or a
% user-supplied Parent axes) can change the scaling; limits stay fixed
% so the view doesn't rescale between frames.
axis(ax, lim*[-1 1 -1 1 -1 1]);
axis(ax, 'equal');
ax.XLimMode = 'manual';  ax.YLimMode = 'manual';  ax.ZLimMode = 'manual';
if showSphere
    axis(axS, 'equal');
    axS.XLimMode = 'manual';  axS.YLimMode = 'manual';  axS.ZLimMode = 'manual';
end

% ---------- Play ----------
if strlength(string(opt.VideoFile)) > 0
    % Write every frame to video
    dt  = median(diff(t));
    fps = min(max(opt.Speed/dt, 1), 120);
    vw  = VideoWriter(char(opt.VideoFile), 'MPEG-4');
    vw.FrameRate = fps;
    open(vw);
    cleanup = onCleanup(@() close(vw));
    for k = 1:N
        if ~isgraphics(hBody), break; end
        drawFrame(k);
        drawnow;
        writeVideo(vw, getframe(fig));
    end
else
    % Real-time playback; skip frames if rendering falls behind
    t0 = t(1);  tEnd = t(end);
    clock0 = tic;
    k = 1;
    while k <= N && isgraphics(hBody)
        drawFrame(k);
        drawnow limitrate;
        simTime = t0 + opt.Speed*toc(clock0);
        if simTime < t(k)                   % ahead of schedule: wait
            pause((t(k) - simTime)/opt.Speed);
            simTime = t(k);
        end
        if simTime >= tEnd && k < N
            k = N;                          % make sure the last frame is shown
        else
            k = max(k+1, find(t <= simTime, 1, 'last') + 1);
        end
    end
    drawnow;                                % flush the final frame
end
end
