function fe = setup_fe(nelx,nely,nu)
% SETUP_FE  Element stiffness, connectivity and sparse assembly indices.

fe.KE = element_stiffness(nu);

nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
fe.edofMat = repmat(edofVec,1,8) + ...
             repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);

fe.iK = reshape(kron(fe.edofMat,ones(8,1))',64*nelx*nely,1);
fe.jK = reshape(kron(fe.edofMat,ones(1,8))',64*nelx*nely,1);
end