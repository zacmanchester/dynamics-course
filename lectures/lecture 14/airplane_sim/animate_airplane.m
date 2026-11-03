function animate_airplane(xhist,xrefhist,thist,p)
%ANIMATE_AIRPLANE Animate the airplane trajectory.
%   The airplane is drawn as a fuselage, a wing, a horizontal tail and a
%   vertical fin. The camera follows the airplane.

    fig = figure('Name','Airplane animation');
    ax = axes(fig);
    hold(ax,'on');
    grid(ax,'on');
    axis(ax,'equal');
    axis(ax,'vis3d');
    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');
    zlabel(ax,'z [m]');
    title(ax,'Airplane animation');
    view(ax,35,25);

    % Grass, runway with dashed centreline, reference path and flown path.
    r = xhist(1:3,:);
    xg = [min(r(1,:))-100, max(r(1,:))+100];
    yg = [min(r(2,:))-100, max(r(2,:))+100];
    patch(ax,xg([1 2 2 1]),yg([1 1 2 2]),zeros(1,4),[0.55 0.75 0.45], ...
          'EdgeColor','none');

    onGround = r(3,:) <= 0;
    if r(3,1) <= 0                                          % flight starts on the runway
        xr = [min(r(1,onGround))-20, max(r(1,onGround))+30];   % runway ends
        wr = 12;                                                % runway width
        patch(ax,xr([1 2 2 1]),0.5*wr*[-1 -1 1 1],0.01*ones(1,4),[0.3 0.3 0.3], ...
              'EdgeColor','none');
        xdash = xr(1):8:xr(2);                                  % 4 m dashes every 8 m
        xdash = [xdash; xdash+4; nan(size(xdash))];
        plot3(ax,xdash(:),zeros(numel(xdash),1),0.02*ones(numel(xdash),1), ...
              'w-','LineWidth',2);
    end
    plot3(ax,xrefhist(1,:),xrefhist(2,:),xrefhist(3,:),'r--');
    traj = animatedline(ax,'LineWidth',1.2);

    % Body geometry in the body frame (x forward, y left, z up), drawn at
    % twice the real size so that the airplane is visible.
    scale = 2;
    Lf = 0.6*p.b;                                        % fuselage length
    fuselage = [0.4*Lf, -0.6*Lf;  0,       0;        0, 0];
    wing     = [0,       0;       0.5*p.b, -0.5*p.b; 0, 0];
    htail    = [-0.6*Lf, -0.6*Lf; 0.2*p.b, -0.2*p.b; 0, 0];
    vfin     = [-0.6*Lf, -0.6*Lf; 0,       0;        0, 0.12*p.b];
    parts  = {scale*fuselage, scale*wing, scale*htail, scale*vfin};
    colors = {'k-','b-','k-','k-'};

    bodyLines = gobjects(4,1);
    for i = 1:4
        bodyLines(i) = plot3(ax,nan(1,2),nan(1,2),nan(1,2),colors{i},'LineWidth',2);
    end

    follow = 20;      % half-size of the viewing box around the airplane [m]
    step   = 5;       % draw every 5th sample

    for k = 1:step:size(xhist,2)
        rk = xhist(1:3,k);
        qk = xhist(4:7,k);
        qk = qk/norm(qk);
        Qk = AirplaneMath.qtoQ(qk);

        for i = 1:4
            pts = rk + Qk*parts{i};
            set(bodyLines(i),'XData',pts(1,:),'YData',pts(2,:),'ZData',pts(3,:));
        end
        addpoints(traj,rk(1),rk(2),rk(3));

        % Follow the airplane with a fixed-size viewing box.
        xlim(ax,rk(1) + [-follow follow]);
        ylim(ax,rk(2) + [-follow follow]);
        zlim(ax,[max(-5,rk(3)-follow), rk(3)+follow]);

        drawnow;

        if k + step <= size(xhist,2)
            pause(max(0,thist(k+step)-thist(k)));
        end
    end
end
