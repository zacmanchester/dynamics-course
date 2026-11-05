function animate_rocket(xhist,uhist,thist,p)
%ANIMATE_ROCKET Animate the rocket landing.
%   The rocket is a thick line along its axis; the flame points along the
%   gimballed thrust direction with a length proportional to the thrust.

    fig = figure('Name','Rocket animation','WindowState','maximized');
    ax = axes(fig);
    hold(ax,'on');
    grid(ax,'on');
    axis(ax,'equal');
    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');
    zlabel(ax,'z [m]');
    title(ax,'Rocket animation');

    r = xhist(1:3,:);
    pad = 2*p.ell;
    xlim(ax,[min(r(1,:))-pad, max(r(1,:))+pad]);
    ylim(ax,[min(r(2,:))-pad, max(r(2,:))+pad]);
    zlim(ax,[-p.ell, max(r(3,:))+pad]);
    view(ax,3);
    rotate3d(ax,'on');                          % drag with the mouse to rotate

    % Landing pad and reference line.
    theta = linspace(0,2*pi,36);
    fill3(ax,p.ell*cos(theta),p.ell*sin(theta),zeros(size(theta)),[0.6 0.6 0.6]);
    plot3(ax,r(1,:),r(2,:),r(3,:),'--');

    body  = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'k-','LineWidth',4);
    flame = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'-','Color',[1 0.5 0],'LineWidth',3);

    bodyPts = [0, 0; 0, 0; -p.ell, p.ell];      % engine to nose, body frame
    step = 5;                                   % draw every 5th sample

    tStart = tic;
    for k = 1:step:size(xhist,2)
        rk = xhist(1:3,k);
        qk = xhist(4:7,k);
        qk = qk/norm(qk);
        Qk = RocketMath.qtoQ(qk);

        pb = rk + Qk*bodyPts;
        set(body,'XData',pb(1,:),'YData',pb(2,:),'ZData',pb(3,:));

        T = uhist(1,k); alpha = uhist(2,k); beta = uhist(3,k);
        dirT = [sin(alpha); sin(beta)*cos(alpha); cos(beta)*cos(alpha)];
        flameLen = p.ell*T/(p.m*p.g);
        pf = rk + Qk*[bodyPts(:,1), bodyPts(:,1) - flameLen*dirT];
        set(flame,'XData',pf(1,:),'YData',pf(2,:),'ZData',pf(3,:));

        drawnow;

        % Pace the animation to wall-clock time.
        pause(max(0,thist(k) - toc(tStart)));
    end
end
