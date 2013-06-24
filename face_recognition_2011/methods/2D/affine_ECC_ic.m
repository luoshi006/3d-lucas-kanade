function fit  = affine_ECC_ic(img, tmplt, p_init, n_iters, verbose, step_size)
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


% Initialisation and pre-computable things
if verbose
    img1 = img; tmplt1 = tmplt;
end
init_a;
tmplt = tmplt - mean(tmplt(:));
n_tmplt = norm(tmplt(:));
tmplt = tmplt/n_tmplt;

% Gradient of template
[tx ty] = gradient(tmplt);

% Jacobian
dW_dp = jacobian_a(w, h);
G = image_jacobian(tx, ty, dW_dp, N_p);

% Hessian and its inverse
C = G' * G;
i_C = inv(C);

% Other precomputable
Gt = G' * tmplt(:);

% Inverse Compositional Algorithm  -------------------------------
for f=1:n_iters
    % Warped image with current parameters
    try
        IWxp = warp_a(img, warp_p, tmplt_pts);
    catch ME
        break;
    end
    IWxp = IWxp-mean(IWxp(:)); % zero-mean image; is useful for brithness change compensation, otherwise you can comment this line
    n_IWxp = norm(IWxp(:));
    IWxp = IWxp/n_IWxp;
    
    % -- Save current fit parameters --
    fit(f).warp_p = warp_p;
    
    % -- Show fitting? --
    if verbose
        IWxp1 = warp_a(img1, warp_p, tmplt_pts);
        error_img = IWxp1 - tmplt1;
        verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
    end
    
    % -- Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    Gw = G' * IWxp(:);
    num = (norm(IWxp(:))^2 - Gw' * i_C * Gw);
    den = (dot(tmplt(:),IWxp(:)) - Gt' * i_C * Gw);
    lambda = num / den;
    
    if den<0 break; end
    
    % Error
    imerror = lambda * IWxp - tmplt;
    Ge = G' * imerror(:);
    
    % Gradient descent parameter updates
    delta_p = i_C * Ge;
    
    % Update warp parmaters
    warp_p = update_step(warp_p, delta_p);
end