function M = build_3d_warp_a(warp_p)
% BUILD_3D_WARP_A - 3D affine warp matrix
%
% [1 -vz vy tx
%  vz 1 -vx ty
% -vy vx 1  tz
%  0  0  0  1]

M = [1         -warp_p(3)  warp_p(2)  warp_p(4); ...
     warp_p(3)  1         -warp_p(1)  warp_p(5); ...
    -warp_p(2)  warp_p(1)  1          warp_p(6); ...
     0          0          0          1];
end