function warp_p = update_step_3d(warp_p, delta_p)
% Compute and apply the update

delta_M = build_3d_warp_a(delta_p);

% Invert compositional warp
delta_M = inv(delta_M);

% Current warp
warp_M = build_3d_warp_a(warp_p);

% Compose
comp_M = warp_M * delta_M;

% Get new parameters - [vx, vy, vz, tx, ty, tz]
warp_p = [comp_M(3,2), comp_M(1,3), comp_M(2,1), ...
          comp_M(1,4), comp_M(2,4), comp_M(3,4)];