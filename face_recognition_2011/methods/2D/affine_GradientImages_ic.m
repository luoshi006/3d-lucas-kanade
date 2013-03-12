function fit = affine_GradientImages_ic(img, tmplt, p_init, n_iters, verbose)
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

% Initialisation and pre-computable things
if verbose
   img1 = img; tmplt1 = tmplt;     
end
init_a;
H  = fspecial('gaussian', [5 5], 2.0);
img = imfilter(img, H, 'replicate');
tmplt = imfilter(tmplt, H, 'replicate');

% Cootes-Taylor feature extraction from template
[nabla_Tx nabla_Ty] = gradient(tmplt);
ab = abs(nabla_Tx + 1i*nabla_Ty); m_ab = ab(:); m_ab = median(m_ab);
ab = ab + m_ab;
tmpltx = nabla_Tx./ab; tmplty = nabla_Ty./ab;

[nabla_Txx nabla_Txy] = gradient(tmpltx);
[nabla_Tyx nabla_Tyy] = gradient(tmplty);
nabla_Tyx = nabla_Txy;

% Jacobian
dW_dp = jacobian_a(w, h);
Gx = image_jacobian(nabla_Txx, nabla_Txy, dW_dp, N_p);
Gy = image_jacobian(nabla_Tyx, nabla_Tyy, dW_dp, N_p);

% Hessian and its inverse
Q = Gx' * Gx + Gy' * Gy;
inv_Q = inv(Q);

% Inverse Compositional Algorithm
for f=1:n_iters
    % Compute warped image and extract Cootes-Taylor features
    try
    IWxp = warp_a(img, warp_p, tmplt_pts);
    catch ME
        break;
    end
    [nabla_Ix nabla_Iy] = gradient(IWxp);
    ab = abs(nabla_Ix + 1i*nabla_Iy); m_ab = ab(:); m_ab = median(m_ab);
    ab = ab + m_ab;
    IWxpx =  nabla_Ix./ab; IWxpy = nabla_Iy./ab;
    
    % Error image
    error_imgx = IWxpx - tmpltx;
    error_imgy = IWxpy - tmplty;
    
    % Save current fit parameters --
    fit(f).warp_p = warp_p;
    
    % Show fitting? --
    if verbose
        IWxp1 = warp_a(img1, warp_p, tmplt_pts);
        error_img = IWxp1 - tmplt1;
        verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
    end
    
    % Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    % Gradient descent parameter updates
    delta_p = inv_Q * (Gx'*error_imgx(:) + Gy'*error_imgy(:));
    
    % Update warp parmaters
    warp_p = update_step(warp_p, delta_p);
end
