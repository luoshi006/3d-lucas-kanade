function wimg = warp_3d_a(img, p, dst)
% WARP_A - Affine warp the image
%   WIMG = WARP_A(IMG, P, DST)
%   Warp image IMG to WIMG. DST are the destination points, i.e. the corners
%   of the template image. P are the affine warp parameters that project
%   DST into IMG.
%
%   P = [p1, p3, p5
%        p2, p4, p6,
%        p7, p8, p9];

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: warp_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<3 error('Not enough input arguments'); end

% Convert affine warp parameters into 4 x 4 warp matrix
% NB affine parameterised as [1 -p1 p2 p3
%                             p4 1 -p5 p6
%                            -p7 p8 1  p9
%                             0  0  0  1]

M = [1   -p(1) p(2)  p(3); ...
     p(4) 1    p(5) -p(6); ...
    -p(7) p(8) 1     p(9); ...
     0    0    0     1];

% Use bilinear filtering to warp image back to template
wimg = polytocuboid(img, dst, M);
