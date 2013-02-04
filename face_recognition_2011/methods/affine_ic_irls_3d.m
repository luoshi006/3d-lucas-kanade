function fitt = affine_ic_irls_3d(img, tmplt, p_init, n_iters, verbose, var)
% AFFINE_IC_IRLS - Affine image alignment using inverse-compositional 
% iteratively reweighted least squares algorithm
%   FIT = AFFINE_IC_IRLS(IMG, TMPLT, P_INIT, N_ITERS, VAR, VERBOSE)
%   Align the template image TMPLT to an example image IMG using an
%   affine warp initialised using P_INIT. Iterate for N_ITERS iterations.
%   VAR.PERC_OUT percent of the pixels are assumed to be outliers.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p3, p5
%             p2, p4, p6];
%
%   This assumes greyscale images and rectangular templates.


% Ralph Gross, Iain Matthews, Simon Baker
% Carnegie Mellon University, Pittsburgh

% Modifications by G. Tzimiropoulos, S. Zafeiriou and M.Pantic 
% To facilitate function call from my_test_affine, the function interface(line 1) has been changed from 
% FIT = AFFINE_IC_IRLS(IMG, TMPLT, P_INIT, N_ITERS, VAR, VERBOSE)
% to
% FIT = AFFINE_IC_IRLS(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE, VAR)
% var.perc_out is manually set (line 34)
% image smoothing (lines 35,36,37) 

if nargin<5 verbose = 0; end
if nargin<4 error('Not enough input arguments'); end

% Common initialisation
[img, warp_p, tmplt_pts, w, h, d, N_p, verb_info] = init_3d_a(tmplt, img, p_init, verbose);

% Pre-computable things ---------------------------------------------------
var.perc_out = 0.4;
% Filter with Gaussian kernel
img = smooth_img(img);
tmplt = smooth_img(tmplt);

% 3) Evaluate gradient of T
[nabla_Tx, nabla_Ty, nabla_Tz] = gradient(tmplt);

% 4) Evaluate Jacobian - constant for affine warps
dW_dp = jacobian_3d_a(w, h, d);

% 5) Compute steepest descent images, VT_dW_dp
VT_dW_dp = sd_images_3d(dW_dp, nabla_Tx, nabla_Ty, nabla_Tz, N_p, w, h, d);

% Algorithm -------------------------
for f=1:n_iters
  % 1) Compute warped image with current parameters
  IWxp = warp_3d_a(img, warp_p, tmplt_pts);

  % 2) Compute error image - NB reversed
  error_img = IWxp - tmplt;
 
  % -- Save current fit parameters --
  fitt(f).warp_p = warp_p;
  fitt(f).rms_error = sqrt(mean(error_img(:) .^2));
  
  % -- Show fitting? --
  if verbose
    verb_plot_3d_a(verb_info, warp_p, tmplt_pts);
  end
  
  % -- Really iteration 1 is the zeroth, ignore final computation --
  if (f == n_iters) break; end

  % Compute robust error funtion
  weight = robust_error_3d(error_img, var.perc_out);

  % 6) Compute robust Hessian
  H     = hessian_weight_3d(VT_dW_dp, weight, N_p, w);
  H_inv = inv(H);

  % 7) Compute weighted steepest descent parameter updates
  sd_delta_p = sd_update_weight_3d(VT_dW_dp, error_img, weight, N_p, w);

  % 8) Compute gradient descent parameter updates
  delta_p = H_inv * sd_delta_p;
  
  % 9) Update warp parmaters
  warp_p = update_step_3d(warp_p, delta_p);
  
end