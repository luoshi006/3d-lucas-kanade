function M = build_3d_warp_2drigid(warp_p)
% BUILD_3D_WARP_2DRIGID - 3D video warp matrix

M = [ warp_p(1) warp_p(3) 0          warp_p(6);
      warp_p(2) warp_p(4) 0          warp_p(7);
      0         0         warp_p(5)  warp_p(8); 
      0         0         0          1         ];

M(1, 1) = M(1, 1) + 1;
M(2, 2) = M(2, 2) + 1;
M(3, 3) = M(3, 3) + 1;