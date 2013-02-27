function fitt  = affine_ECC_ic_3d(img, tmplt, p_init, n_iters, verbose, smoothing)
% affine_ECC_ic - Affine image alignment using the enhanced correlation coefficient [1] and the
% inverse-compositional algorithm of Baker-Matthews [2]
%
%   FIT = affine_ECC_ic(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
%   Align the template image TMPLT to an example image IMG using an
%   affine warp initialised using P_INIT. Iterate for N_ITERS iterations.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p3, p5
%             p2, p4, p6];
%
%   This assumes greyscale images and rectangular templates.
%
% References:
% [1] G.D. Evangelidis and  E.Z. Psarakis, "Parametric Image Alignment using Enhanced Correlation Coefficient Maximization", IEEE TPAMI, Vol. 30, No. 10, pp. 1858-1865, 2008
% [2] S. Baker and I. Matthews. "Equivalence and effciency of image alignment algorithms", Proc. CVPR 2001, pp 1090-1097,
%
% Implemented using functions and code provided by 
% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% http://www.ri.cmu.edu/research_project_detail.html?project_id=515&menu_id=261
% AND
% Georgios Evangelidis, Visual and Social Media Lab (VSM), Fraunhofer IAIS, St. Augustin, Germany
% http://xanthippi.ceid.upatras.gr/people/evangelidis/

% Written by G. Tzimiropoulos, S. Zafeiriou and M. Pantic 
% Intelligent Behaviour Understanding Group (IBUG), Department of Computing, Imperial College London
% $ Version: 1.0, 03/01/2012 $


if nargin<5 verbose = 0; end
if nargin<4 error('Not enough input arguments'); end

% Common initialisation
[img, warp_p, tmplt_pts, w, h, d, N_p, verb_info] = init_3d_a(tmplt, img, p_init, verbose);

% Filter with Gaussian kernel
if (smoothing)
    img = smooth_img(img);
    tmplt = smooth_img(tmplt);
end

tmplt   = tmplt - mean(tmplt(:));
n_tmplt = norm(tmplt(:));
tmplt   = tmplt / n_tmplt;

% 3) Evaluate gradient of T
[nabla_Tx, nabla_Ty, nabla_Tz] = gradient(tmplt);

% 4) Evaluate Jacobian - constant for affine warps
dW_dp = jacobian_3d_a(w, h, d);
G     = image_jacobian_3d(nabla_Tx, nabla_Ty, nabla_Tz, dW_dp, N_p);

% Hessian and its inverse
C   = G' * G;
i_C = inv(C);

% Other precomputable
Gt = G' * tmplt(:);

% Inverse Compositional Algorithm  -------------------------------
for f=1:n_iters
    % Warped image with current parameters
    try
        IWxp = warp_3d_a(img, warp_p, tmplt_pts);
    catch ME
        break;
    end
    
    % zero-mean image is useful for brightness change compensation
    IWxp   = IWxp - mean(IWxp(:)); % otherwise you can comment this line
    n_IWxp = norm(IWxp(:));
    IWxp   = IWxp / n_IWxp;
    
    % -- Save current fit parameters --
    fitt(f).warp_p = warp_p;
    
    % -- Show fitting? --
    if verbose
        verb_plot_3d_a(verb_info, warp_p, tmplt_pts);
    end
    
    % -- Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    Gw     = G' * IWxp(:);
    num    = (norm(IWxp(:)) ^ 2 - Gw' * i_C * Gw);
    den    = (dot(tmplt(:), IWxp(:)) - Gt' * i_C * Gw);
    lambda = num / den;
    
    if den < 0 
        fprintf('The denominator is diverging: %f \n', den);
        break; 
    end
    
    % Error
    imerror = lambda * tmplt - IWxp;
    Ge      = G' * imerror(:);
    
    % Gradient descent parameter updates
    delta_p = -i_C * Ge;
    
    % Update warp parmaters
    warp_p = update_step_3d(warp_p, delta_p);
end