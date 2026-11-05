function [pen,n] = ground(pw,p)
%GROUND Tire compression and ground normal at the world point pw.
%   pen > 0 means the point is pushed into the ground by pen meters, and n
%   is the unit normal of the ground there (pointing out of the ground).
%
%   The ground is the flat plane z = 0. A demo may switch on a loop in the
%   x-z plane: a quarter circle of radius R that starts flat at x0 and turns
%   vertical, a vertical wall of height Lv, a second quarter circle that
%   curves back over the car, and a straight inverted segment of length Lt
%   at the top, after which the track ends.

    pen = -pw(3);
    n   = [0; 0; 1];

    if p.loop.on
        ztop = 2*p.loop.R + p.loop.Lv;
        if pw(3) > p.loop.R && pw(3) <= p.loop.R + p.loop.Lv && pw(1) > p.loop.x0
            % Vertical wall at x = x0 + R.
            pen = pw(1) - (p.loop.x0 + p.loop.R);
            n   = [-1; 0; 0];
        elseif pw(3) > p.loop.R + p.loop.Lv && pw(1) <= p.loop.x0 && pw(1) >= p.loop.x0 - p.loop.Lt
            % Straight inverted segment at z = ztop.
            pen = pw(3) - ztop;
            n   = [0; 0; -1];
        else
            % Lower or upper circle, depending on the height.
            zc = p.loop.R + p.loop.Lv*(pw(3) > p.loop.R + p.loop.Lv);
            d  = [pw(1) - p.loop.x0; pw(3) - zc];       % from the circle center
            theta = atan2(d(1),-d(2));                  % 0 at the bottom, pi at the top
            if theta >= 0
                rho = norm(d);
                pen = rho - p.loop.R;                   % beyond the track
                n   = -[d(1); 0; d(2)]/rho;             % toward the center
            end
        end
    end
end
