function wimg = warp_3d_a(img, p, dst)
% WARP_A - Affine warp the image
%   WIMG = WARP_A(IMG, P, DST)
%   Warp image IMG to WIMG. DST are the destination points, i.e. the corners
%   of the template image. P are the affine warp parameters that project
%   DST into IMG.
%
%   P = [p1, p3, p5   = [vx, vy, vz
%        p2, p4, p6]    tx, ty, tz]

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: warp_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<3 error('Not enough input arguments'); end

% Convert affine warp parameters into 4 x 4 warp matrix
% NB affine parameterised as [1 -vz vy tx
%                             vz 1 -vx ty
%                            -vy vx 1  tz
%                             0  0  0  1]

M = build_3d_warp_a(p);

% Use bilinear filtering to warp image back to template
wimg = polytocuboid(img, dst, M);
