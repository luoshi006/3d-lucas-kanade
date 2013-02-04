function M = build_3d_warp_a(warp_p)
% BUILD_3D_WARP_A - 3D affine warp matrix

warp_p = reshape(warp_p, 3, 4);

M = [warp_p; 0 0 0 1];
M(1, 1) = M(1, 1) + 1;
M(2, 2) = M(2, 2) + 1;
M(3, 3) = M(3, 3) + 1;
end