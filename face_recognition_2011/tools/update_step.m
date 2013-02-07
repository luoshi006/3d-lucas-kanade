function warp_p = update_step(warp_p, delta_p)
% Compute and apply the update

delta_p = reshape(delta_p, 2, 3);

% Convert affine notation into usual Matrix form - NB transposed
delta_M = [delta_p; 0 0 1];
delta_M(1,1) = delta_M(1,1) + 1;
delta_M(2,2) = delta_M(2,2) + 1;

% Invert compositional warp
delta_M = inv(delta_M);

% Current warp
warp_M = [warp_p; 0 0 1];
warp_M(1,1) = warp_M(1,1) + 1;
warp_M(2,2) = warp_M(2,2) + 1;

% Compose
comp_M = warp_M * delta_M;
warp_p = comp_M(1:2,:);
warp_p(1,1) = warp_p(1,1) - 1;
warp_p(2,2) = warp_p(2,2) - 1;