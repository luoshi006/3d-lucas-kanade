function verb_plot_3d_a(verb, warp_p, tmplt_pts, error_img)
% VERB_PLOT_A - Verbose fitting plot
%   VERB_PLOT_A(V, WARP_P, TMPLT_PTS, ERROR_IMG)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: verb_plot_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<2 error('Not enough input arguments'); end

% Current parameters
warp_M = build_3d_warp_a(warp_p);
warp_pts =  warp_M * [tmplt_pts; ones(1, size(tmplt_pts,2))];
set(verb.lh, 'Xdata', [warp_pts(1,:) warp_pts(1,1)], ...
             'Ydata', [warp_pts(2,:) warp_pts(2,1)], ...
             'Zdata', [warp_pts(3,:) warp_pts(3,1)]);
drawnow
axis off