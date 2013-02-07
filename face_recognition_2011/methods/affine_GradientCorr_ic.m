function fitt = affine_GradientCorr_ic(img, tmplt, p_init, n_iters, verbose)
% affine_GradientCorr_ic - Affine image alignment using the proposed maximization of gradient correlation [1] 
%
%   FIT = affine_GradientCorr_ic(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
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
% [1] "Robust and effiecient face alignment by gradient correlation maximization", submitted to ICCV 2011
%
% Based on the implementation of
% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: affine_ic.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $
% http://www.ri.cmu.edu/research_project_detail.html?project_id=515&menu_id=261
% AND
% Georgios Evangelidis, Visual and Social Media Lab (VSM), Fraunhofer IAIS, St. Augustin, Germany
% http://www.mathworks.co.uk/matlabcentral/fileexchange/27253-ecc-image-alignment-algorithm

if nargin<5 verbose = 0; end
if nargin<4 error('Not enough input arguments'); end


% Initialisation and pre-computable things 
if verbose
   img1 = img; 
   tmplt1 = tmplt;     
end
init_a;
H     = fspecial('gaussian', [5 5], 2.0);
img   = imfilter(img, H, 'replicate');
tmplt = imfilter(tmplt, H, 'replicate');
% Gradients of image - will warp these later
%[vx1 vy1] = myEdge(img, 5);
% vabs = abs(vx1 + 1i*vy1);
% vx1 = vx1./vabs; vy1 = vy1./vabs;
% vx1(isnan(vx1))=0;
% vy1(isnan(vy1))=0;

% Gradient of template
[tx, ty] = myEdge(tmplt, 5);
ang = angle(tx + 1i * ty);
tx = cos(ang); 
ty = sin(ang);

[txx, txy] = myEdge(tx, 5);
[tyx, tyy] = myEdge(ty, 5);
tyx = txy;

fx = -sin(ang); 
fx = fx(:); 
fx(isnan(fx)) = 0; 
fx = repmat(fx, 1, N_p); 

fy = cos(ang);
fy = fy(:); 
fy(isnan(fy)) = 0;  
fy = repmat(fy, 1, N_p);

% Jacobian 
dW_dp = jacobian_a(w, h);
Gx = image_jacobian(txx, txy, dW_dp, N_p);
Gy = image_jacobian(tyx, tyy, dW_dp, N_p);

G = fx .* Gx + fy .* Gy;
Gxx = fx .* Gx; 
Gyy = fy .* Gy;

% Hessian and its inverse
Q = Gxx' * Gxx + Gyy' * Gyy;

% Other precomputable
v = numel(tmplt);

% Inverse Compositional Algorithm  -------------------------------
for f=1:n_iters
    % Warped image with current parameters (don't warp first iteration)
    if (exist('warp_p', 'var'))
        IWxp = warp_a(img, warp_p, tmplt_pts);  
    end
    
    [vx, vy] = myEdge(IWxp, 5);
    ang = angle(vx + 1i * vy);
    vx = cos(ang); 
    vy = sin(ang);

    
    % -- Save current fit parameters --
    fitt(f).warp_p = warp_p;
    
    % -- Show fitting? --
    if verbose
        IWxp1 = warp_a(img1, warp_p, tmplt_pts);
        error_img = IWxp1 - tmplt1;
        verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
    end
    
    % -- Really iteration 1 is the zeroth, ignore final computation --
    if (f == n_iters) break; end
    
    u_bold = G' * (tx(:) .* vy(:) - ty(:) .* vx(:));
    u      = vx(:)' * tx(:) + vy(:)' * ty(:);
    
    % lambda 
    lambda = v / u;
    if u < 0 
        break; 
    end
    
    % Error 
    imerror = lambda * u_bold;
    
    % Gradient descent parameter updates
    %delta_p = inv_Q * imerror;
    delta_p =  Q \ imerror;

    % Update warp parmaters
    warp_p = update_step(warp_p, delta_p);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bx, by] = myEdge(I, par)

if par == 1
    gx = [-1, 0, 1] / 2;
elseif par == 2
    gx = [1/12, -2/3, 0, 2/3, -1/12];
elseif par == 3
    gx = [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60];
elseif par == 4
    sigma = 1;
    gx = Gradx_oG(max(1, floor(5 * sigma)), sigma);
elseif par == 5
    gx = -dxmask;
elseif par == 6
    gx = -dxxmask;
elseif par == 7
    gx = -dxymask;
end
    
bx = crop2(conv2(extend2(I, length(gx), length(gx)), gx, 'same'), length(gx), length(gx));     
by = crop2(conv2(extend2(I, length(gx), length(gx)), gx', 'same'), length(gx), length(gx));   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = dxmask
result = [0 0 0  0 0 0 0 0 0; 
          0 0 0  0 0 0 0 0 0;  
          0 0 0  0 0 0 0 0 0;  
          0 0 0  0 0 0 0 0 0; 
          0 0 0 -1 0 1 0 0 0 
          0 0 0  0 0 0 0 0 0; 
          0 0 0  0 0 0 0 0 0; 
          0 0 0  0 0 0 0 0 0; 
          0 0 0  0 0 0 0 0 0] / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = crop2(data, ny, nx)
[ysize, xsize] = size(data);
result = data((ny + 1:ysize - ny), (nx + 1:xsize - nx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = extend2(data, ny, nx)
[ysize, xsize] = size(data);
newxsize = xsize + 2 * nx;
newysize = ysize + 2 * ny;

result = zeros(newysize, newxsize);
result((ny + 1:ysize + ny), (nx + 1:xsize + nx)) = data;

for x = 1:nx
    result(:, x) = result(:, nx + 1); 
end

for x = (xsize + nx + 1):newxsize
    result(:, x) = result(:, xsize + nx); 
end

for y = 1:ny
    result(y, :) = result(ny + 1, :); 
end

for y = (ysize + ny + 1):newysize
    result(y, :) = result(ysize + ny, :); 
end
