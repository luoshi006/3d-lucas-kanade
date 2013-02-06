function weight = robust_error_3d(error_img, perc_out)
% ROBUST_ERROR - Robust error function
%  WEIGHT = ROBUST_ERROR (ERROR_IMG, PERC_OUT)
%
% Compute robust error function of the error image ERROR_IMG
% assuming that PERC_OUT percent of the pixels are outliers

% Ralph Gross
% Carnegie Mellon University, Pittsburgh

[~, ind]       = sort(abs(error_img(:)));
selInd         = ind(1:round((1 - perc_out) * length(ind)));
weight         = zeros(size(error_img));
weight(selInd) = 1;


