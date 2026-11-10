
function U = FE_analysis(nelx,nely,x,penal,E0,Emin,...
                         KE,F,freedofs)

ndof = 2*(nelx+1)*(nely+1);

% Element connectivity
nele = nelx*nely;
edofMat = zeros(nele,8);

e = 0;

for j = 1:nely
    for i = 1:nelx

        e = e+1;

        n1 = (j-1)*(nelx+1)+i;
        n2 = n1+1;
        n4 = j*(nelx+1)+i;
        n3 = n4+1;

        edofMat(e,:) = ...
            [2*n1-1,2*n1,...
             2*n2-1,2*n2,...
             2*n3-1,2*n3,...
             2*n4-1,2*n4];
    end
end

% Assemble global stiffness matrix
iK = reshape(kron(edofMat,ones(8,1))',[],1);
jK = reshape(kron(edofMat,ones(1,8))',[],1);

sK = zeros(64*nele,1);

for e = 1:nele

    rho = x(e);

    Ee = Emin + rho^penal*(E0-Emin);

    sK(64*(e-1)+1:64*e) = ...
        reshape(Ee*KE,64,1);

end

K = sparse(iK,jK,sK,ndof,ndof);

% Solve equilibrium
U = zeros(ndof,1);

U(freedofs) = K(freedofs,freedofs)\F(freedofs);

end
