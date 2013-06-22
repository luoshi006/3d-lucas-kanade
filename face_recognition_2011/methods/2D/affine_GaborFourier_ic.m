function fit = affine_GaborFourier_ic(img, tmplt, p_init, n_iters, verbose, S)
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

if nargin<5 verbose = 0; end
if nargin<4 error('Not enough input arguments'); end

% Initialisation and pre-computable things
if verbose
   img1 = img; tmplt1 = tmplt;     
end
init_a;
SS = repmat(S, 1, N_p);

% Gradient of template
[nabla_Tx, nabla_Ty] = gradient(tmplt);

% Jacobian 
dW_dp = jacobian_a(w, h);

% Fourier Steepest descent images, VT_dW_dp
VT_dW_dp = sd_images(dW_dp, nabla_Tx, nabla_Ty, N_p, h, w);
FT_VT_dW_dp = zeros(size(VT_dW_dp));
len = size(tmplt,2);
for ii = 1:N_p
    FT_VT_dW_dp(:,(ii-1)*len+1:ii*len ) = fftshift(fft2(VT_dW_dp(:,(ii-1)*len+1:ii*len )));
end

% Hessian and inverse
H = my_hessian((SS.^(1/2)).*FT_VT_dW_dp, N_p, w);
H_inv = inv(H);

% Gabor Fourier Inverse Compositional Algorithm 
for f=1:n_iters
    % Warped image with current parameters
    try
    IWxp = warp_a(img, warp_p, tmplt_pts);
    catch ME
        break;
    end
    
    % Error image 
    error_img = IWxp - tmplt;
    FT_error_img = S.*fftshift(fft2(error_img));
    
    % Save current fit parameters 
    fit(f).warp_p = warp_p;
    
    % -- Show fitting? --
    if verbose
        IWxp1 = warp_a(img1, warp_p, tmplt_pts);
        error_img = IWxp1 - tmplt1;
        verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
    end
    
    % -- Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    % Steepest descent parameter updates
    sd_delta_p = sd_update(FT_VT_dW_dp, FT_error_img, N_p, w);

    % Gradient descent parameter updates
    delta_p = H_inv * sd_delta_p;
    delta_p = real(delta_p);
    
    % 9) Update warp parmaters
    warp_p = update_step(warp_p, delta_p);
    warp_p = real(warp_p);
end
