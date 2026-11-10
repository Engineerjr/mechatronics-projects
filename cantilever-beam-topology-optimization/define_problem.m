function prob = define_problem(nelx,nely)
% DEFINE_PROBLEM  Loads and boundary conditions for the cantilever beam.
%   Left edge fully clamped; unit downward point load at mid-height of the
%   right edge. Edit this file to change the load case.
%
%   DOF numbering: node n has DOFs (2n-1 = x, 2n = y); nodes are numbered
%   column-wise from the top-left corner.

prob.nelx = nelx;
prob.nely = nely;
prob.ndof = 2*(nely+1)*(nelx+1);

% Load: downward unit force at mid-height of right edge
loadNode = nelx*(nely+1) + round(nely/2) + 1;
prob.F = sparse(2*loadNode,1,-1,prob.ndof,1);

% Supports: clamp every DOF on the left edge
prob.fixeddofs = 1:2*(nely+1);
prob.freedofs  = setdiff(1:prob.ndof,prob.fixeddofs);
end