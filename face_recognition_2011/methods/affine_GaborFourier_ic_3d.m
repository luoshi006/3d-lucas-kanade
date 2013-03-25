function fitt = affine_GaborFourier_ic_3d(img, tmplt, p_init, n_iters, verbose, smoothing, varargin)
% affine_GaborFourier_ic - Affine image alignment using the Gabor-Fourier framework [1] and the
% inverse-compositional algorithm of Baker-Matthews [2]
%
%   FIT = affine_GaborFourier_ic(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
%   Align the template image TMPLT to an example image IMG using an
%   affine warp initialised using P_INIT. Matrix S is described in [1]. Iterate for N_ITERS iterations.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p3, p5
%             p2, p4, p6];
%
%   This assumes greyscale images and rectangular templates.
%
% References:
% [1] A.B. Ashraf, S. Lucey and T. Chen, "Fast Image Alignment in the Fourier Domain", Proc. CVPR 2010, pp 2480-2487 
% [2] S. Baker and I. Matthews. "Equivalence and effciency of image alignment algorithms", Proc. CVPR 2001, pp 1090-1097,
%
% Implemented using functions and code provided by 
% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% http://www.ri.cmu.edu/research_project_detail.html?project_id=515&menu_id=261
% AND
% Peter Kovesi's Gabor Filters
% http://www.csse.uwa.edu.au/~pk/

% Written by G. Tzimiropoulos, S. Zafeiriou and M. Pantic 
% Intelligent Behaviour Understanding Group (IBUG), Department of Computing, Imperial College London
% $ Version: 1.0, 03/01/2012 $

% Common initialisation
[img, warp_p, tmplt_pts, w, h, d, N_p, verb_info] = init_3d_a(tmplt, img, p_init, verbose);
S = varargin{1};
SS = repmat(S, 1, N_p);

% Gradient of template
[nabla_Tx, nabla_Ty, nabla_Tz] = gradient(tmplt);

% Jacobian 
dW_dp = jacobian_3d_a(w, h, d);

% Fourier Steepest descent images, VT_dW_dp
VT_dW_dp = sd_images_3d(dW_dp, nabla_Tx, nabla_Ty, nabla_Tz, N_p, w, h, d);
FT_VT_dW_dp = zeros(size(VT_dW_dp));
for p = 1:N_p
    index = ((p - 1) * w) + 1:((p - 1) * w) + w;
    FT_VT_dW_dp(:, index, :) = fftshift(fftn(VT_dW_dp(:, index, :)));
end

% Hessian and inverse
H = hessian_3d(sqrt(SS) .* FT_VT_dW_dp, N_p, w, true);
H_inv = inv(H);

% Gabor Fourier Inverse Compositional Algorithm 
for f=1:n_iters
    % Warped image with current parameters
    IWxp = warp_3d_a(img, warp_p, tmplt_pts);
    
    % Error image 
    error_img = IWxp - tmplt;
    FT_error_img = S .* fftshift(fftn(error_img));
    
    % Save current fit parameters 
    fitt(f).warp_p = warp_p;
    
    % -- Show fitting? --
    if verbose
    	verb_plot_3d_a(verb_info, warp_p, tmplt_pts);
    end
    
    % -- Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    % Steepest descent parameter updates
    sd_delta_p = sd_update_3d(FT_VT_dW_dp, FT_error_img, N_p, w);

    % Gradient descent parameter updates
    delta_p = H_inv * sd_delta_p;
    delta_p = real(delta_p);
    
    % 9) Update warp parmaters
    warp_p = update_step_3d(warp_p, delta_p);
    warp_p = real(warp_p);
end
