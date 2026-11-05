function animate_flyingcar(xhist,uhist,thist,p)
%ANIMATE_FLYINGCAR Animate the flying car.
%   The car is a red box with four wheels. Orange flames show the thrusters
%   (length proportional to thrust) and a blue line below the car shows
%   the suction force. The loop and the obstacle cars are drawn when the
%   demo defines them (p.loop and p.cars).

    fig = figure('Name','Flying car animation','WindowState','maximized');
    ax = axes(fig);
    hold(ax,'on');
    grid(ax,'on');
    axis(ax,'equal');
    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');
    zlabel(ax,'z [m]');
    title(ax,'Flying car animation');

    r = xhist(1:3,:);
    pad = 4;
    xlim(ax,[min(r(1,:))-pad, max(r(1,:))+pad]);
    ylim(ax,[min(r(2,:))-pad, max(r(2,:))+pad]);
    zlim(ax,[-1, max(r(3,:))+pad]);
    view(ax,-35,15);
    rotate3d(ax,'on');                          % drag with the mouse to rotate

    % Flat ground.
    xg = [min(r(1,:))-pad, max(r(1,:))+pad];
    yg = [-pad, pad];
    [Xg,Yg] = meshgrid(xg,yg);
    surf(ax,Xg,Yg,zeros(2),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none');

    % Loop: a ribbon along the track.
    if p.loop.on
        [xl,zl] = loop_profile(p);
        surf(ax,[xl; xl],[-2*ones(size(xl)); 2*ones(size(xl))],[zl; zl], ...
             'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.3 0.3 0.3]);
    end

    % Obstacle cars: grey boxes.
    if isfield(p,'cars')
        for i = 1:numel(p.cars.x)
            draw_box(ax,[p.cars.x(i); 0; p.cars.H(i)/2],p.cars.L,p.cars.W,p.cars.H(i),[0.4 0.4 0.4]);
        end
    end

    % Body box, wheels and flames.
    L = p.a + 0.6;                              % half length of the body
    boxPts = [ L,  L, -L, -L,  L,  L, -L, -L; ...
               p.b, -p.b, -p.b,  p.b,  p.b, -p.b, -p.b,  p.b; ...
              -0.2, -0.2, -0.2, -0.2, 0.4, 0.4, 0.4, 0.4];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    body = patch(ax,'Vertices',boxPts','Faces',faces, ...
                 'FaceColor',[0.85 0.1 0.1],'FaceAlpha',0.9);

    rwheel = 0.35;
    theta = linspace(0,2*pi,24);
    circle = rwheel*[cos(theta); zeros(size(theta)); sin(theta)];
    wheels = gobjects(4,1);
    flames = gobjects(4,1);
    for i = 1:4
        wheels(i) = plot3(ax,nan(size(theta)),nan(size(theta)),nan(size(theta)), ...
                          'k-','LineWidth',2);
        flames(i) = plot3(ax,nan(1,2),nan(1,2),nan(1,2), ...
                          '-','Color',[1 0.5 0],'LineWidth',3);
    end
    flame5  = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'-','Color',[1 0.5 0],'LineWidth',3);
    suction = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'-','Color',[0 0.4 1],'LineWidth',3);

    step = 5;                                   % draw every 5th sample

    tStart = tic;
    for k = 1:step:size(xhist,2)
        rk = xhist(1:3,k);
        qk = xhist(4:7,k);
        qk = qk/norm(qk);
        Qk = FlyingCarMath.qtoQ(qk);

        pb = rk + Qk*boxPts;
        set(body,'Vertices',pb');

        for i = 1:4
            cw = p.rw(:,i) + [0; 0; rwheel];                 % wheel center
            pw = rk + Qk*(cw + circle);
            set(wheels(i),'XData',pw(1,:),'YData',pw(2,:),'ZData',pw(3,:));

            flameLen = 2*uhist(i,k)/p.Tmax;
            pf = rk + Qk*[p.rw(:,i), p.rw(:,i) - [0; 0; flameLen]];
            set(flames(i),'XData',pf(1,:),'YData',pf(2,:),'ZData',pf(3,:));
        end

        flameLen = 2*uhist(5,k)/p.T5max;
        pf = rk + Qk*[-L, -L - flameLen; 0, 0; 0, 0];
        set(flame5,'XData',pf(1,:),'YData',pf(2,:),'ZData',pf(3,:));

        suctionLen = 1*uhist(6,k)/p.Smax;
        ps = rk + Qk*[0, 0; 0, 0; -0.2, -0.2 - suctionLen];
        set(suction,'XData',ps(1,:),'YData',ps(2,:),'ZData',ps(3,:));

        drawnow;

        % Pace the animation to wall-clock time.
        pause(max(0,thist(k) - toc(tStart)));
    end
end

function draw_box(ax,c,L,W,H,color)
    % Axis-aligned box centered at c with length L, width W and height H.
    pts = c + [L/2*[1 1 -1 -1 1 1 -1 -1]; ...
               W/2*[1 -1 -1 1 1 -1 -1 1]; ...
               H/2*[-1 -1 -1 -1 1 1 1 1]];
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    patch(ax,'Vertices',pts','Faces',faces,'FaceColor',color,'FaceAlpha',0.8);
end
