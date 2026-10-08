function body = makeTBody(varargin)
%MAKETBODY Geometry and inertia of a uniform T-shaped rigid body.
%
%   body = makeTBody() returns a T made of two uniform rectangular boxes:
%   a crossbar and a stem. All quantities are expressed in the body frame,
%   with the origin at the center of mass:
%       x : thickness direction (normal to the plane of the T)
%       y : along the crossbar
%       z : along the stem (stem points toward -z)
%   By symmetry these are principal axes, so body.J is diagonal.
%
%   Name-value options (consistent units):
%     'BarLength'   crossbar length            (default: 1.0)
%     'StemLength'  stem length                (default: 0.8)
%     'Width'       in-plane width of the bars (default: 0.15)
%     'Thickness'   out-of-plane thickness     (default: 0.15)
%     'Density'     mass density               (default: 1)
%
%   Output struct fields:
%     J          3x3 inertia matrix about the center of mass
%     mass       total mass
%     Vertices   16x3 patch vertices (body frame, COM at origin)
%     Faces      12x4 patch faces
%     FacePart   12x1 part index for each face (1 = crossbar, 2 = stem)
%     tip        3x1 body-frame position of the end of the stem
%     com        center of mass in the original geometry frame
%     dims       struct of the dimensions used
%
%   Example:
%     body = makeTBody('StemLength', 1.2);
%     patch('Vertices',body.Vertices, 'Faces',body.Faces, 'FaceColor','c');
%     axis equal; view(3)
%
%   See also ANIMATERIGIDBODY.

p = inputParser;
pos = @(x) isnumeric(x) && isscalar(x) && x > 0;
p.addParameter('BarLength', 1.0, pos);
p.addParameter('StemLength', 0.8, pos);
p.addParameter('Width', 0.15, pos);
p.addParameter('Thickness', 0.15, pos);
p.addParameter('Density', 1, pos);
p.parse(varargin{:});
d = p.Results;

Lb = d.BarLength;  Ls = d.StemLength;  w = d.Width;  th = d.Thickness;

% ---------- Parts: [center; full side lengths] in geometry frame ----------
% Crossbar centered at the origin; stem hangs below it along -z.
cBar  = [0, 0, 0];               sBar  = [th, Lb, w];
cStem = [0, 0, -w/2 - Ls/2];     sStem = [th, w, Ls];

mBar  = d.Density*prod(sBar);
mStem = d.Density*prod(sStem);
mTot  = mBar + mStem;
com   = (mBar*cBar + mStem*cStem)/mTot;

% ---------- Inertia about the COM (box + parallel axis theorem) ----------
J = boxInertia(mBar, sBar, cBar - com) + boxInertia(mStem, sStem, cStem - com);

% ---------- Patch geometry (shifted so the COM is at the origin) ----------
[vB, fB] = boxPatch(cBar - com, sBar);
[vS, fS] = boxPatch(cStem - com, sStem);

body.J        = J;
body.mass     = mTot;
body.Vertices = [vB; vS];
body.Faces    = [fB; fS + size(vB,1)];
body.FacePart = [ones(size(fB,1),1); 2*ones(size(fS,1),1)];
body.tip      = ([0, 0, -w/2 - Ls] - com).';
body.com      = com.';
body.dims     = d;
end

% =====================================================================
function J = boxInertia(m, s, r)
% Inertia of a uniform box (side lengths s) about a point offset by -r
% from its center, i.e. the box center sits at position r.
Jc = (m/12)*diag([s(2)^2 + s(3)^2, s(1)^2 + s(3)^2, s(1)^2 + s(2)^2]);
r  = r(:);
J  = Jc + m*((r.'*r)*eye(3) - (r*r.'));
end

function [V, F] = boxPatch(c, s)
% Vertices and quad faces of a box with center c and side lengths s.
u = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
V = c + 0.5*u.*s;
F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
end
