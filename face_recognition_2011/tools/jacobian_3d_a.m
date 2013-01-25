function dW_dp = jacobian_3d_a(nx, ny, nz);
% JACOBIAN_A - Compute Jacobian for affine warp
%   DW_DP = JACOBIAN_A(WIDTH, HEIGHT, DEPTH)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: jacobian_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% TODO: How do i arrange these?
jac_x = kron((0:nx - 1), ones(ny, 1));
jac_y = kron((0:ny - 1)', ones(1, nx));
jac_z = kron((0:nz - 1)', ones(1, nx));
jac_zero = zeros(ny, nx);
jac_one = ones(ny, nx);

% [ 0  z -y 1 0 0 
%  -z  0  x 0 1 0
%   y -x  0 0 0 1 ]
dW_dp = [
         jac_zero, jac_z, -jac_y, jac_one, jac_zero, jac_zero;
        -jac_z, jac_zero, jac_x, jac_zero, jac_one, jac_zero;
         jac_y, -jac_x, jac_zero, jac_zero, jac_zero, jac_one;
        ];
