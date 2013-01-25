function VI_dW_dp = sd_images_3d(dW_dp, nabla_Ix, nabla_Iy, nabla_Iz, N_p, h, w, d)
% SD_IMAGES - Compute steepest descent images
%   VI_DW_DP = SD_IMAGES(DW_DP, NABLA_IX, NABLA_IY, N_P, H, W)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: sd_images.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<6 error('Not enough input arguments'); end
% TODO: fix this
for p=1:N_p		
	Tx = nabla_Ix .* dW_dp(1:h,((p-1)*w)+1:((p-1)*w)+w);
	Ty = nabla_Iy .* dW_dp(h+1:end,((p-1)*w)+1:((p-1)*w)+w);
	VI_dW_dp(:,((p-1)*w)+1:((p-1)*w)+w) = Tx + Ty;
end
