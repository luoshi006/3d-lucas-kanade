function wimg = polytocuboid(img, dst, M)
% POLYTOCUBOID - Warp contents of a polyhedron to a cuboid
%   W = POLYTOCUBOID(IMG, DST, M, F_TYPE)
%
%   Destination corners, DST = [x1, x2, x3, x4; y1, y2, y3, y4];
%   These are the (assumed rectangular) corner points of the template
%   image and must be integers.
%
%   The affine matrix, M transforms DST into the image IMG 
%   [u v 1] (image) = M * [x y 1] (DST - template);
%
%   Matlab style indices, i.e. start from 1.

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: quadtobox.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<3 error('Invalid args'); end

% Dimensions of destination image - integers, assume cuboid
maxv = max(dst');

% Get all points in destination to sample
[xg, yg, zg] = meshgrid(1:maxv(1), 1:maxv(2), 1:maxv(3));
xy = [xg(:)'; ...
      yg(:)'; ...
      zg(:)'];
xy = [xy; ones(1, size(xy, 2))];

% Transform into source
uv = M * xy;

% Remove homogeneous
uv = uv(1:3,:)';

% Sample
wimg = interp3(img, uv(:,1), uv(:,2), uv(:,3));
wimg = reshape(wimg, [maxv(1) maxv(2) maxv(3)]);

% Check for NaN background pixels - replace them with a background of 0
idx = find(isnan(wimg));
if ~isempty(idx)
	wimg(idx) = 0;
end
