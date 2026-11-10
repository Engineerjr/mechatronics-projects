function [H,Hs] = build_filter(nelx,nely,rmin)
% BUILD_FILTER  Sparse linear (cone) filter matrix and its row sums.

r  = ceil(rmin)-1;
nmax = nelx*nely*(2*r+1)^2;
iH = zeros(nmax,1);
jH = zeros(nmax,1);
sH = zeros(nmax,1);
k  = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely + j1;
        for i2 = max(i1-r,1):min(i1+r,nelx)
            for j2 = max(j1-r,1):min(j1+r,nely)
                e2 = (i2-1)*nely + j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin - sqrt((i1-i2)^2 + (j1-j2)^2));
            end
        end
    end
end
H  = sparse(iH(1:k),jH(1:k),sH(1:k),nelx*nely,nelx*nely);
Hs = sum(H,2);
end