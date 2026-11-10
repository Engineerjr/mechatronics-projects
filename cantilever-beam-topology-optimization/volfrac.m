nelx    = 120;       % Number of elements along X
nely    = 40;        % Number of elements along Y
volfrac = 0.40;      % Maximum material volume fraction

% Initial density field: every element starts at volfrac
x = volfrac * ones(nely, nelx);

% Quantities derived from volfrac
nElem      = nelx * nely;              % total number of elements
targetVol  = volfrac * nElem;          % allowed material volume (in element units)
currentVol = sum(x(:));                % current material volume
currentFrac = currentVol / nElem;      % current volume fraction

fprintf('Elements          : %d x %d = %d\n', nelx, nely, nElem);
fprintf('Target volume     : %.1f elements\n', targetVol);
fprintf('Current volume    : %.1f elements\n', currentVol);
fprintf('Current fraction  : %.3f (limit %.2f)\n', currentFrac, volfrac);

% Plot the density field (white = void, black = solid)
figure;
colormap(gray);
imagesc(1 - x);                        % 1-x so high density appears dark
axis([0 1]);
axis equal tight off;
title(sprintf('Initial design: %d x %d elements, volfrac = %.2f', nelx, nely, volfrac));
colorbar;