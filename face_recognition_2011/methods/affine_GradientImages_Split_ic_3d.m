function fitt = affine_GradientImages_Split_ic_3d(img, tmplt, p_init, n_iters, verbose, smoothing, varargin)
% affine_GradientImages_ic - Affine image alignment using the features proposed by Cootes-Taylor [1] and the
% inverse-compositional algorithm of Baker-Matthews [2]
%
%   FIT = affine_GradientImages_ic(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
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
% [1] T.F.Cootes and C.J.Taylor. "On Representing Edge Structure for Model Matching", Proc. CVPR 2001, Volume 1, pp.1114-1119
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

% Cootes-Taylor feature extraction from template
[g1x, g1y, g1z, g1sz] = median_adjusted_gradient(tmplt);

[nabla_Txx, nabla_Txy, nabla_Txz] = gradient(g1x);
[nabla_Tyx, nabla_Tyy, nabla_Tyz] = gradient(g1y);
[nabla_Tzx, nabla_Tzy, nabla_Tzz] = gradient(g1z);
[nabla_Tszx, nabla_Tszy, nabla_Tszz] = gradient(g1sz);
% Approximates fully differentiable gradients
nabla_Tyx = nabla_Txy;
nabla_Tyz = nabla_Tzy;
nabla_Txz = nabla_Tzx;

% Jacobian
dW_dp = jacobian_3d_a(w, h, d);
Gx = image_jacobian_3d(nabla_Txx, nabla_Txy, nabla_Txz, dW_dp, N_p);
Gy = image_jacobian_3d(nabla_Tyx, nabla_Tyy, nabla_Tyz, dW_dp, N_p);
Gz = image_jacobian_3d(nabla_Tzx, nabla_Tzy, nabla_Tzz, dW_dp, N_p);
Gsz = image_jacobian_3d(nabla_Tszx, nabla_Tszy, nabla_Tszz, dW_dp, N_p);

% Hessian and its inverse
Q     = Gx' * Gx + Gy' * Gy + Gz' * Gz + Gsz' * Gsz;
inv_Q = inv(Q);

% Inverse Compositional Algorithm
for f=1:n_iters
    % Compute warped image and extract Cootes-Taylor features
    try
        IWxp = warp_3d_a(img, warp_p, tmplt_pts);
    catch ME
        break;
    end
    [g2x, g2y, g2z, g2sz] = median_adjusted_gradient(IWxp);
    
    % Error image
    error_imgx = g2x - g1x;
    error_imgy = g2y - g1y;
    error_imgz = g2z - g1z;
    error_imgsz = g2sz - g1sz;
    
    % Save current fit parameters --
    fitt(f).warp_p = warp_p;
    
    % Show fitting? --
    if verbose
        verb_plot_3d_a(verb_info, warp_p, tmplt_pts);
    end
    
    % Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    % Gradient descent parameter updates
    delta_p = inv_Q * (Gx' * error_imgx(:) + Gy' * error_imgy(:) + Gz' * error_imgz(:) + Gsz' * error_imgsz(:));
    
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