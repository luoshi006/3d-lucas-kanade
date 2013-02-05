function sd = image_jacobian_3d(gx, gy, gz, dW_dp, N_p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%G = IMAGE_JACOBIAN(GX, GY, J, NOP)
% This function computes the jacobian G of warped image wrt parameters. 
% This matrix depends on the gradient of warped image, as 
% well as of the jacobian J of the warp transform wrt parameters. 
% For a detailed definition of matrix G, see Evangelidis & Psarakis paper.
%
% Input variables:
% GX:           the warped image gradient in x (horizontal) deirection,
% GY:           the warped image gradient in y (horizontal) deirection,
% DW_DP:            the jacobian matrix J of warp transform wrt parameters,
% N_P:          the number of parameters.
%
% Output:
% G:            The jacobian matrix G.
%--------------------------------------
% $ Ver: 1.0.0, 1/3/2010,  released by Georgios D. Evangelidis, Fraunhofer IAIS.
% For any comment, please contact georgios.evangelidis@iais.fraunhofer.de
% or evagelid@ceid.upatras.gr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h, w, d] = size(gx);

if nargin<4
    error('Not enough input arguments');
end

gx = repmat(gx, 1, N_p);
gy = repmat(gy, 1, N_p);
gz = repmat(gz, 1, N_p);

Tx = gx .* dW_dp(1:h, :, :);
Ty = gy .* dW_dp(h+1:h+h, :, :);
Tz = gz .* dW_dp(h+h+1:end, :, :);
G = Tx + Ty + Tz;

sd = zeros(w * h * d, N_p);
for p=1:N_p
    cols = ((p - 1) * w) + 1:((p - 1) * w) + w;
    sd(:, p) = reshape(G(:, cols, :), [w * h * d, 1]);
end