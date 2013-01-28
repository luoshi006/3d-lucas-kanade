function VI_dW_dp = sd_images_3d(dW_dp, nabla_Ix, nabla_Iy, nabla_Iz, N_p, w, h, d)
% SD_IMAGES - Compute steepest descent images
%   VI_DW_DP = SD_IMAGES(DW_DP, NABLA_IX, NABLA_IY, N_P, H, W)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: sd_images.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<6 error('Not enough input arguments'); end

VI_dW_dp = zeros(h, w * N_p, d);
for p=1:N_p		
    cols = ((p-1)*w)+1:((p-1)*w)+w;
	Tx = nabla_Ix .* dW_dp(1:h, cols, :);
	Ty = nabla_Iy .* dW_dp(h+1:h+h, cols, :);
    Tz = nabla_Iz .* dW_dp(h+h+1:end, cols, :);
	VI_dW_dp(:, cols, :) = Tx + Ty + Tz;
end
