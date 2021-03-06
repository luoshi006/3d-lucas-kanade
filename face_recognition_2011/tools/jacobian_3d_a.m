function dW_dp = jacobian_3d_a(nx, ny, nz)
% JACOBIAN_A - Compute Jacobian for affine warp
%   DW_DP = JACOBIAN_A(WIDTH, HEIGHT, DEPTH)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: jacobian_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

[jac_x, jac_y, jac_z] = meshgrid(0:nx-1, 0:ny-1, 0:nz-1);
jac_zero = zeros(ny, nx, nz);
jac_one = ones(ny, nx, nz);

dW_dp = [
         jac_x,    jac_zero, jac_zero, jac_y,    jac_zero, jac_zero, jac_z,    jac_zero, jac_zero, jac_one,  jac_zero, jac_zero;
         jac_zero, jac_x,    jac_zero, jac_zero, jac_y,    jac_zero, jac_zero, jac_z,    jac_zero, jac_zero, jac_one,  jac_zero
         jac_zero, jac_zero, jac_x,    jac_zero, jac_zero, jac_y,    jac_zero, jac_zero, jac_z,    jac_zero, jac_zero, jac_one;
        ];
