function sd_delta_p = sd_update_weight_3d(VI_dW_dp, error_img, weight, N_p, w)
% SD_UPDATE_WEIGHT - Compute weighted steepest descent parameter updates
%   SD_DELTA_P = SD_UPDATE_WEIGHT(VI_DW_DP, ERROR_IMG, WEIGHT, N_P, W)

% Iain Matthews, Ralph Gross,
% Carnegie Mellon University, Pittsburgh

if nargin<5 error('Not enough input arguments'); end

sd_delta_p = zeros(N_p, 1);
for p=1:N_p
  h1 = VI_dW_dp(:, ((p - 1) * w) + 1:((p - 1) * w) + w, :);
  sd_delta_p(p) = sum(sum(sum(weight .* (h1 .* error_img))));
end
