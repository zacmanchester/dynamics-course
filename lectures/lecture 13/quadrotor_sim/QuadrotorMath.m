classdef QuadrotorMath
%QUADROTORMATH Small collection of quaternion and rigid-body helpers.

    methods (Static)
        function S = hat(v)
            S = [ 0,    -v(3),  v(2); ...
                  v(3),  0,    -v(1); ...
                 -v(2),  v(1),  0];
        end

        function M = Lquat(q)
            s = q(1);
            v = q(2:4);
            M = [s, -v'; ...
                 v, s*eye(3) + QuadrotorMath.hat(v)];
        end

        function H = Hmat()
            H = [zeros(1,3); eye(3)];
        end

        function Q = qtoQ(q)
            T = diag([1,-1,-1,-1]);
            H = QuadrotorMath.Hmat();
            L = QuadrotorMath.Lquat(q);
            Q = H'*T*L*T*L*H;
        end

        function G = Gquat(q)
            G = QuadrotorMath.Lquat(q)*QuadrotorMath.Hmat();
        end

        function M = E(q)
            M = blkdiag(eye(3),QuadrotorMath.Gquat(q),eye(6));
        end

        function q = rptoq(phi)
            q = [1; phi]/sqrt(1 + phi'*phi);
        end

        function phi = qtorp(q)
            phi = q(2:4)/q(1);
        end
    end
end
