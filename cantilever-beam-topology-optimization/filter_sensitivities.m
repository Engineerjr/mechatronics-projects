function [dc,dv] = filter_sensitivities(dc,dv,x,H,Hs,ft)
% FILTER_SENSITIVITIES  Apply sensitivity (ft=1) or density (ft=2) filtering.

switch ft
    case 1   % sensitivity filter
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    case 2   % density filter (chain rule)
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    otherwise
        error('ft must be 1 (sensitivity) or 2 (density).');
end
end