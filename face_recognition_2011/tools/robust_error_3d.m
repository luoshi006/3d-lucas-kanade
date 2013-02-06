function weight = robust_error_3d(error_img, k)
% ROBUST_ERROR - Robust error function
%  WEIGHT = ROBUST_ERROR (ERROR_IMG, PERC_OUT)
%
% Compute a Huber M-estimator

% Ralph Gross
% Carnegie Mellon University, Pittsburgh

weight = zeros(size(error_img));
error_img = abs(error_img);
weight(error_img <= k) = 1;
weight(error_img > k) = k / error_img(error_img > k);


