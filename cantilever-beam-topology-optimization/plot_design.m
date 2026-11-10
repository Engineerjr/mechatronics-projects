function plot_design(xPhys)
% PLOT_DESIGN  Display the current density field (black = solid).

colormap(gray);
imagesc(1-xPhys);
caxis([0 1]);
axis equal; axis off;
drawnow;
end