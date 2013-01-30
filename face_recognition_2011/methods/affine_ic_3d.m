function fitt = affine_ic_3d(img, tmplt, p_init, n_iters, verbose, ~)
% AFFINE_IC - Affine image alignment using inverse-compositional algorithm
%   FIT = AFFINE_IC(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
%   Align the template image TMPLT to an example image IMG using an
%   affine warp initialised using P_INIT. Iterate for N_ITERS iterations.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p3, p5
%             p2, p4, p6];
%
%   This assumes greyscale images and rectangular templates.
%
%   c.f. Baker-Matthews

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: affine_ic.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

% Image smoothing (lines 27,28,29) added by G. Tzimiropoulos, S. Zafeiriou and M. Pantic

if nargin<5 
    verbose = 0; 
end
if nargin<4 
    error('Not enough input arguments'); 
end

% Common initialisation
[img, warp_p, tmplt_pts, w, h, d, verb_info] = init_3d_a(tmplt, img, p_init, verbose);

% Pre-computable things ---------------------------------------------------
% TODO: do this in 3D
%H  = fspecial('gaussian', [5 5], 2.0);
%img = imfilter(img, H, 'replicate'); 
%tmplt = imfilter(tmplt, H, 'replicate');

% 3) Evaluate gradient of T
[nabla_Tx, nabla_Ty, nabla_Tz] = gradient(tmplt);

% 4) Evaluate Jacobian - constant for affine warps
dW_dp = jacobian_3d_a(w, h, d);

% 5) Compute steepest descent images, VT_dW_dp
VT_dW_dp = sd_images_3d(dW_dp, nabla_Tx, nabla_Ty, nabla_Tz, N_p, w, h, d);
	
% 6) Compute Hessian and inverse
H = hessian_3d(VT_dW_dp, N_p, w);
H_inv = inv(H);

% Baker-Matthews, Inverse Compositional Algorithm -------------------------

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
        verb_plot_3d_a(verb_info, warp_p, tmplt_pts, error_img);
	end
	
	% -- Really iteration 1 is the zeroth, ignore final computation --
	if (f == n_iters) 
        break;
    end

	% 7) Compute steepest descent parameter updates
	sd_delta_p = sd_update_3d(VT_dW_dp, error_img, N_p, w);

	% 8) Compute gradient descent parameter updates
	delta_p = H_inv * sd_delta_p;
	
	% 9) Update warp parmaters
	warp_p = update_step_3d(warp_p, delta_p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function warp_p = update_step_3d(warp_p, delta_p)
% Compute and apply the update

delta_M = build_3d_warp_a(delta_p);

% Invert compositional warp
delta_M = inv(delta_M);

% Current warp
warp_M = build_3d_warp_a(warp_p);

% Compose
comp_M = warp_M * delta_M;

% Get new parameters - [vx, vy, vz, tx, ty, tz]
warp_p = [comp_M(3,2), comp_M(1,3), comp_M(2,1), ...
          comp_M(1,4), comp_M(2,4), comp_M(3,4)];
