function fitt = affine_GradientCorr_Euclidean_ic(img, tmplt, p_init, n_iters, verbose)
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
init_a;

% Filter with Gaussian kernel
H     = fspecial('gaussian', [5 5], 2.0);
img   = imfilter(img, H, 'replicate');
tmplt = imfilter(tmplt, H, 'replicate');

% 3) Evaluate gradient of T
[g1x, g1y] = gradient(tmplt);

% Calculate df(g1,x[0]) / dg1,x[0]
g1x2 = g1x .^ 2;
g1y2 = g1y .^ 2;

g1_norm = abs(g1x + 1i*g1y);
df_g1_denom = g1_norm .^ 3;

m_ab = median(g1_norm(:));
g1_norm = g1_norm + m_ab;
df_g1_denom = df_g1_denom + m_ab;

dF_g1x = g1y2 ./ df_g1_denom;
dF_g1y = g1x2 ./ df_g1_denom;

tilde_g1x = g1x ./ g1_norm;
tilde_g1y = g1y ./ g1_norm;
G1 = tilde_g1x + tilde_g1y;
G1 = G1(:);

% Calculate dg1,x[0]/dp
[g1xx, g1xy] = gradient(g1x);
[g1yx, g1yy] = gradient(g1y);
g1yx = g1xy;

% [g1,xx, g1,xy] * dW_dp
dW_dp = jacobian_a(w, h);
dg1x_dp = image_jacobian(g1xx, g1xy, dW_dp, N_p);
dg1y_dp = image_jacobian(g1yx, g1yy, dW_dp, N_p);

dF_g1x = repmat(dF_g1x(:), 1, N_p);
dF_g1y = repmat(dF_g1y(:), 1, N_p);
Jx = dF_g1x .* dg1x_dp;
Jy = dF_g1y .* dg1y_dp;

J = Jx + Jy;

% Hessian and its inverse
H = Jx' * Jx + Jy' * Jy;
invH = inv(H);

% Inverse Compositional Algorithm  -------------------------------
for f=1:n_iters
    % Warped image with current parameters
    IWxp = warp_a(img, warp_p, tmplt_pts);
    
    % -- Save current fit parameters --
    fitt(f).warp_p = warp_p;
    
    [g2x, g2y] = gradient(IWxp);
    g2_norm = abs(g2x + 1i*g2y);
    g2_norm = g2_norm + median(g2_norm(:));
    tilde_g2x = g2x ./ g2_norm;
    tilde_g2y = g2y ./ g2_norm;
    G2 = tilde_g2x + tilde_g2y;
    G2 = G2(:);
    
    % -- Show fitting? --
    if verbose
        verb_plot_a(verb_info, warp_p, tmplt_pts, IWxp - tmplt);
    end
    
    % -- Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    Gw      = J' * G1;
%   wPgw    = Gw' * invH * Gw;
    Gr      = J' * G2;
%   rPgr    = Gr' * invH * Gr;
%   rPgw    = Gr' * invH * Gw;
%     
%   lambda1 = sqrt(wPgw / rPgr);
%   lambda2 = (rPgw - dot(G2, G1)) / rPgr;
%   lambda  = max(lambda1, lambda2);
    
    num    = norm(G1)^2 - Gw' * invH * Gw;
    den    = dot(G2, G1) - Gr' * invH * Gw;
    lambda = num / den;

    if den < 0 
        fprintf('The denominator is diverging: %f', den);
        break; 
    end
    
    % Error
    imerrorx = lambda * tilde_g2x(:) - tilde_g1x(:);
    imerrory = lambda * tilde_g2y(:) - tilde_g1y(:);
    Ge      = Jx' * imerrorx + Jy' * imerrory;
    
    % Gradient descent parameter updates
    delta_p = invH * Ge;
    
    % Update warp parmaters
    warp_p = update_step(warp_p, delta_p);
end