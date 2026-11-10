function [c,dc,dv] = compute_sensitivities(xPhys,U,fe,mat)
% COMPUTE_SENSITIVITIES  Compliance and its derivative, plus volume gradient.

[nely,nelx] = size(xPhys);
ce = reshape(sum((U(fe.edofMat)*fe.KE).*U(fe.edofMat),2),nely,nelx);

c  = sum(sum((mat.Emin + xPhys.^mat.penal*(mat.E0-mat.Emin)).*ce));
dc = -mat.penal*(mat.E0-mat.Emin)*xPhys.^(mat.penal-1).*ce;
dv = ones(nely,nelx);
end