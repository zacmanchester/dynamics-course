function body = makeBoxBody(varargin)
%MAKEBOXBODY Geometry and inertia of a uniform rectangular box.
%
%   body = makeBoxBody() returns a uniform box centered at its center of
%   mass, with edges along the body x, y, z axes (its principal axes).
%   With distinct side lengths the three principal moments are distinct:
%   the longest side is the minimum-inertia axis, the shortest side the
%   maximum-inertia axis.
%
%   Name-value options (consistent units):
%     'Lx', 'Ly', 'Lz'  side lengths              (default: 1.2, 0.8, 0.4)
%     'Density'         mass density              (default: 1)
%
%   Output struct fields (same layout as makeTBody, for animateRigidBody):
%     J          3x3 inertia matrix about the center of mass (diagonal)
%     mass       total mass
%     Vertices   8x3 patch vertices (body frame)
%     Faces      6x4 patch faces
%     FacePart   6x1 part index per face (1 = +/-x, 2 = +/-y, 3 = +/-z)
%     tip        3x1 center of the +x face (traced by the animation trail)
%     dims       struct of the dimensions used
%
%   See also MAKETBODY, ANIMATERIGIDBODY.

p = inputParser;
pos = @(x) isnumeric(x) && isscalar(x) && x > 0;
p.addParameter('Lx', 1.2, pos);
p.addParameter('Ly', 0.8, pos);
p.addParameter('Lz', 0.4, pos);
p.addParameter('Density', 1, pos);
p.parse(varargin{:});
d = p.Results;

s = [d.Lx, d.Ly, d.Lz];
m = d.Density*prod(s);

body.J    = (m/12)*diag([s(2)^2 + s(3)^2, s(1)^2 + s(3)^2, s(1)^2 + s(2)^2]);
body.mass = m;

u = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
body.Vertices = 0.5*u.*s;
% Faces ordered: -x, +x, -y, +y, -z, +z
body.Faces    = [1 4 8 5; 2 3 7 6; 1 2 6 5; 4 3 7 8; 1 2 3 4; 5 6 7 8];
body.FacePart = [1; 1; 2; 2; 3; 3];
body.tip      = [s(1)/2; 0; 0];
body.dims     = d;
end
