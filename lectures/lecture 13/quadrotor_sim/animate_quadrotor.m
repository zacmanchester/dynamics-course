function animate_quadrotor(xhist,thist,armLength)
%ANIMATE_QUADROTOR Animate the quadrotor trajectory.
%   Blue circles at the four arm endpoints represent the propellers.

    fig = figure('Name','Quadrotor animation');
    ax = axes(fig);
    hold(ax,'on');
    grid(ax,'on');
    axis(ax,'equal');
    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');
    zlabel(ax,'z [m]');
    title(ax,'Quadrotor animation');

    r = xhist(1:3,:);
    propRadius = 0.22*armLength;
    pad = max(0.5,2*armLength + propRadius);
    xlim(ax,[min(r(1,:))-pad, max(r(1,:))+pad]);
    ylim(ax,[min(r(2,:))-pad, max(r(2,:))+pad]);
    zlim(ax,[min(0,min(r(3,:))-pad), max(r(3,:))+pad]);
    view(ax,3);

    plot3(ax,r(1,:),r(2,:),r(3,:),'--');
    arm1 = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'-','LineWidth',2);
    arm2 = plot3(ax,nan(1,2),nan(1,2),nan(1,2),'-','LineWidth',2);

    bodyPts1 = armLength*[1, -1; 0, 0; 0, 0];
    bodyPts2 = armLength*[0, 0; 1, -1; 0, 0];
    rotorCentersBody = armLength*[1, -1, 0, 0; ...
                                  0,  0, 1,-1; ...
                                  0,  0, 0, 0];

    % Propeller circles lie in the quadrotor body x-y plane.
    theta = linspace(0,2*pi,36);
    circleBody = propRadius*[cos(theta); sin(theta); zeros(size(theta))];

    props = gobjects(4,1);
    for i = 1:4
        props(i) = plot3(ax,nan(size(theta)),nan(size(theta)),nan(size(theta)), ...
                         'b-','LineWidth',1.5);
    end

    for k = 1:size(xhist,2)
        rk = xhist(1:3,k);
        qk = xhist(4:7,k);
        qk = qk/norm(qk);
        Qk = QuadrotorMath.qtoQ(qk);

        p1 = rk + Qk*bodyPts1;
        p2 = rk + Qk*bodyPts2;

        set(arm1,'XData',p1(1,:),'YData',p1(2,:),'ZData',p1(3,:));
        set(arm2,'XData',p2(1,:),'YData',p2(2,:),'ZData',p2(3,:));

        for i = 1:4
            rotorCircle = rk + Qk*(rotorCentersBody(:,i) + circleBody);
            set(props(i), ...
                'XData',rotorCircle(1,:), ...
                'YData',rotorCircle(2,:), ...
                'ZData',rotorCircle(3,:));
        end

        drawnow;

        if k < size(xhist,2)
            pause(max(0,thist(k+1)-thist(k)));
        end
    end
end
