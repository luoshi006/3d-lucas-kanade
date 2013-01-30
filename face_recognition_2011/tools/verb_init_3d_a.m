function verb = verb_init_3d_a(img, tmplt, tmplt_pts, warp_p)
% VERB_INIT_A - Initialise verbose plot
%   V = VERB_INIT_A(IMG, TMPLT, TMPLT_PTS, WARP_P)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: verb_init_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% Modified by by G. Tzimiropoulos, S. Zafeiriou and M. Pantic 

% TODO: fix this

if nargin<4 error('Not enough input arguments'); end

% Init figure
clf;
set(gcf,'DoubleBuffer','on');
colormap(gray(256));

% Input image
subplot(2,2,1);
patch(isosurface(img,0.99));
axis('image');
title('Image');
axis off


warp_M = build_3d_warp_a(warp_p);
warp_pts = warp_M * [tmplt_pts; ones(1, size(tmplt_pts, 2))];

hold on
verb.lh = plot3([warp_pts(1,:) warp_pts(1,1)], [warp_pts(2,:) warp_pts(2,1)], [warp_pts(3,:) warp_pts(3,1)], 'linewidth', 2);
hold off

% Template image
subplot(2,2,2);
patch(isosurface(tmplt, 0.99));
axis('image');
title('Template');
axis off

% Error image
subplot(2,2,4);
verb.ih_error = patch(isosurface(tmplt, 0.99));
axis('image');
title('Error');
axis off
drawnow;
