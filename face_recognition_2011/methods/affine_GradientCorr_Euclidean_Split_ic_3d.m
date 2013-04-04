function fitt = affine_GradientCorr_Euclidean_Split_ic_3d(img, tmplt, p_init, n_iters, verbose, smoothing, varargin)
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

% 3) Evaluate gradient of T
% Need both tilde_g1 and g1
[g1x, g1y, g1z] = gradient(tmplt);
[tilde_g1x, tilde_g1y, tilde_g1z, tilde_g1sz] = median_adjusted_gradient(tmplt);

% Calculate df(g1,x[0]) / dg1,x[0]
g1x2 = g1x .^ 2;
g1y2 = g1y .^ 2;
g1z2 = g1z .^ 2;

df_g1xy_denom = sqrt(g1x2 + g1y2) .^ 3;
df_g1z_denom = sqrt(g1x2 + g1y2 + g1z2) .^ 3;
df_g1sz_denom1 = g1x2 + g1y2 + g1z2;
df_g1sz_denom2 = g1y2 + g1z2;

% Prevent division by zero by adding the median
% Take the median ignoring the dead voxels
df_g1xy_denom = median_adjusted(df_g1xy_denom);
df_g1z_denom = median_adjusted(df_g1z_denom);
df_g1sz_denom1 = median_adjusted(df_g1sz_denom1);
df_g1sz_denom2 = median_adjusted(df_g1sz_denom2);

dF_g1x = g1y2 ./ df_g1xy_denom;
dF_g1y = g1x2 ./ df_g1xy_denom;
dF_g1z = (g1x2 + g1y2) ./ df_g1z_denom;
dF_g1sz = -((g1z .* ((g1y2 + g1z2) ./ df_g1sz_denom1) .^ (3/2)) ./ df_g1sz_denom2);

% Calculate dg1,x[0]/dp
[g1xx, g1xy, g1xz] = gradient(g1x);
[g1yx, g1yy, g1yz] = gradient(g1y);
[g1zx, g1zy, g1zz] = gradient(g1z);
g1yx = g1xy;
g1yz = g1zy;
g1xz = g1zx;

% [g1,xx, g1,xy, g1,xz] * dW_dp
dW_dp = jacobian_3d_a(w, h, d);
dg1x_dp = image_jacobian_3d(g1xx, g1xy, g1xz, dW_dp, N_p);
dg1y_dp = image_jacobian_3d(g1yx, g1yy, g1yz, dW_dp, N_p);
dg1z_dp = image_jacobian_3d(g1zx, g1zy, g1zz, dW_dp, N_p);

dF_g1x = repmat(dF_g1x(:), 1, N_p);
dF_g1y = repmat(dF_g1y(:), 1, N_p);
dF_g1z = repmat(dF_g1z(:), 1, N_p);
dF_g1sz = repmat(dF_g1sz(:), 1, N_p);
Jx = dF_g1x .* dg1x_dp;
Jy = dF_g1y .* dg1y_dp;
Jz = dF_g1z .* dg1z_dp;
Jsz = dF_g1sz .* dg1z_dp;

% Hessian and its inverse
H = Jx' * Jx + Jy' * Jy + Jz' * Jz + Jsz' * Jsz;
invH = inv(H);

% Inverse Compositional Algorithm  -------------------------------
for f=1:n_iters
    % Warped image with current parameters
    IWxp = warp_3d_a(img, warp_p, tmplt_pts);
    
    % -- Save current fit parameters --
    fitt(f).warp_p = warp_p;
    
    % -- Show fitting? --
    if verbose
        verb_plot_3d_a(verb_info, warp_p, tmplt_pts);
    end
    
    % -- Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    [tilde_g2x, tilde_g2y, tilde_g2z, tilde_g2sz] = median_adjusted_gradient(IWxp);
    
    % Error image
    error_imgx = tilde_g2x - tilde_g1x;
    error_imgy = tilde_g2y - tilde_g1y;
    error_imgz = tilde_g2z - tilde_g1z;
    error_imgsz = tilde_g2sz - tilde_g1sz;
    
    % Save current fit parameters --
    fitt(f).warp_p = warp_p;
    
    % Show fitting? --
    if verbose
        verb_plot_3d_a(verb_info, warp_p, tmplt_pts);
    end
    
    % Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    % Gradient descent parameter updates
    delta_p = invH * (Jx' * error_imgx(:) + Jy' * error_imgy(:) + Jz' * error_imgz(:) + Jsz' * error_imgsz(:));
    
    % Update warp parmaters
    warp_p = update_step_3d(warp_p, delta_p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gx, gy, gz, sgz] = median_adjusted_gradient(img)
[nabla_Tx, nabla_Ty, nabla_Tz] = gradient(img);

xyz = sqrt(nabla_Tx.^2 + nabla_Ty.^2 + nabla_Tz.^2);
xyz = median_adjusted(xyz);

xy = sqrt(nabla_Tx.^2 + nabla_Ty.^2);
xy = median_adjusted(xy);

gx  = nabla_Tx ./ xy; 
gy  = nabla_Ty ./ xy;
gz  = nabla_Tz ./ xyz;
sgz = sqrt(1 - gz .^ 2);