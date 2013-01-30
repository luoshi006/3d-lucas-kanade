function [img, warp_p, tmplt_pts, w, h, d, verb_info] = init_3d_a(tmplt, img, p_init, verbose)
% init_a.m
% Common initialisation things for all affine algorithms

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: init_a.m,v 1.1.1.1 2113/18/21 13:17:36 iainm Exp $

% Need to process image data, must be real
if ~isa(img, 'double')
	img = double(img);
end

% 3D Affine warp (skew matrix, translation)
N_p = 6;
if numel(p_init) ~= N_p
	error('Number of warp parameters incorrect');
end

% Initial warp parameters
warp_p = p_init;

% Template size
w = size(tmplt, 1);
h = size(tmplt, 2);
d = size(tmplt, 3);

% Template verticies, cube           
tmplt_pts = [
             1 1 1;
             1 h 1;
             w h 1;
             w 1 1;
             1 1 d;
             1 h d;
             w h d;
             w 1 d;
            ]';

% Verbose display of fitting?
if verbose
	verb_info = verb_init_3d_a(img, tmplt, tmplt_pts, warp_p);
end
