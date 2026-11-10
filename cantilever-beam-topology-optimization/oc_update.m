function [xnew,xPhys] = oc_update(x,dc,dv,H,Hs,volfrac,ft)
% OC_UPDATE  Optimality criteria update with bisection on the multiplier.

[nely,nelx] = size(x);
xPhys = x;
l1 = 0; l2 = 1e9; move = 0.2;

while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if ft == 1
        xPhys = xnew;
    else
        xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nelx*nely
        l1 = lmid;
    else
        l2 = lmid;
    end
end
end