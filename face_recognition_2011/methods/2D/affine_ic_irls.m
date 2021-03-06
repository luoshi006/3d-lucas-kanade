function fit = affine_ic_irls(img, tmplt, p_init, n_iters, verbose, var)
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
init_a;

% Pre-computable things ---------------------------------------------------
var.perc_out = 0.4;
H  = fspecial('gaussian', [5 5], 2.0); 
img = imfilter(img, H, 'replicate');
tmplt = imfilter(tmplt, H, 'replicate');

% 3) Evaluate gradient of T
[nabla_Tx nabla_Ty] = gradient(tmplt);

% 4) Evaluate Jacobian - constant for affine warps
dW_dp = jacobian_a(w, h);

% 5) Compute steepest descent images, VT_dW_dp
VT_dW_dp = sd_images(dW_dp, nabla_Tx, nabla_Ty, N_p, h, w);

% Algorithm -------------------------
for f=1:n_iters
  % 1) Compute warped image with current parameters
  IWxp = warp_a(img, warp_p, tmplt_pts);

  % 2) Compute error image - NB reversed
  error_img = IWxp - tmplt;
 
  % -- Save current fit parameters --
  fit(f).warp_p = warp_p;
  fit(f).rms_error = sqrt(mean(error_img(:) .^2));
  
  % -- Show fitting? --
  if verbose
    %disp(['Inverse-Compositional [',num2str(f-1),']: RMS = ',num2str(fit(f).rms_error)]);
    verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
  end
  
  % -- Really iteration 1 is the zeroth, ignore final computation --
  if (f == n_iters) break; end

  % Compute robust error funtion
  weight = robust_error (error_img, var.perc_out);

  % 6) Compute robust Hessian
  H     = hessian_weight(VT_dW_dp, weight, N_p, w);
  H_inv = inv(H);

  % 7) Compute weighted steepest descent parameter updates
  sd_delta_p = sd_update_weight(VT_dW_dp, error_img, weight, N_p, w);

  % 8) Compute gradient descent parameter updates
  delta_p = H_inv * sd_delta_p;
  
  % 9) Update warp parmaters
  warp_p = update_step(warp_p, delta_p);
  
end
