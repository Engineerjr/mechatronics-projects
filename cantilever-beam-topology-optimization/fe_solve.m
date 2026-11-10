function U = fe_solve(xPhys,fe,prob,mat)
% FE_SOLVE  Assemble global stiffness matrix and solve K*U = F.

n  = prob.nelx*prob.nely;
sK = reshape(fe.KE(:)*(mat.Emin + xPhys(:)'.^mat.penal*(mat.E0-mat.Emin)),64*n,1);
K  = sparse(fe.iK,fe.jK,sK);
K  = (K+K')/2;

U = zeros(prob.ndof,1);
U(prob.freedofs) = K(prob.freedofs,prob.freedofs)\prob.F(prob.freedofs);
end