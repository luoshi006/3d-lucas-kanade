function results = test_affine_3d(tdata, pt_offsets, alg_list, n_iters, n_freq_tests, spatial_sigma, max_spatial_error, verbose, smoothing)
% my_test_affine - Test affine algorithms: a short version of the
% test_affine function written by Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
%
% See also: main_Yale, main_AR, test_affine
%
% tdata has three fields:
%   tdata.img1      unit8 greyscale EXAMPLE IMAGE
%   tdata.img2      unit8 greyscale image used to generate the
%                   TEMPLATE IMAGE
%   tdata.tmplt     [x_start y_start x_end y_end] template rectangle corners
%                   matlab coordinates (1, 1) is top-left, start < end.
%
% pt_offsets is a N x 6 matrix of (x1, y1, x2, y2, x3, y3) deltas for each
% of the three affine control points. Zero mean, unit variance,
% e.g. pt_offsets = randn(N, 6);
%
% alg_list is a cellstr list of algorithms to run, e.g.:
% alg_list = {'affine_ic' 'affine_ECC_ic' 'affine_ic_irls' 'affine_GaborFourier_ic' 'affine_GradientImages_ic' 'affine_GradientCorr_ic'};

% Written by G. Tzimiropoulos, S. Zafeiriou and M. Pantic 
% Intelligent Behaviour Understanding Group (IBUG), Department of Computing, Imperial College London
% $ Version: 1.0, 03/01/2012 $

if nargin<8 error('Not enough input arguments'); end

minX = tdata.tmplt(1);
minY = tdata.tmplt(2);
minZ = tdata.tmplt(3);
maxX = tdata.tmplt(4);
maxY = tdata.tmplt(5);
maxZ = tdata.tmplt(6);

% Target affine warp control points - tetrahedron
target_affine = [minX minY minZ;                                                           % Bottom triangle
                 maxX minY minZ;
                 minX + ((maxX - minX) / 2) - 0.5 maxY minZ;
                 minX + ((maxX - minX) / 2) - 0.5 minY + ((maxY - minY) / 2) - 0.5 maxZ]'; % Tip of tetrahedron

% Template image dimensions
template_nx = maxX - minX;
template_ny = maxY - minY;
template_nz = maxZ - minZ;

% Template corner points (unperturbed, rectangle)
template_pts = [1           1           1;           % Bottom Rectangle
				1           template_ny 1;
				template_nx template_ny 1;
				template_nx 1           1
                1           1           template_nz; % Top Rectangle
				1           template_ny template_nz;
				template_nx template_ny template_nz;
				template_nx 1           template_nz]';

% Template affine warp control points
template_affine = [1               1               1;             % Bottom Triangle
				   template_nx     1               1;
				   template_nx / 2 template_ny     1;
                   template_nx / 2 template_ny / 2 template_nz]'; % Tip of Tetrahedron Triangle

% Initial warp parameters. Unperturbed translation
p_init = zeros(3, 4);
p_init(:, 4) = [minX - 1; minY - 1; minZ - 1];

% Translate by 0.5 pixels to avoid identity warp. Warping introduces a little
% smoothing and this avoids the case where the first iteration for a forwards
% algorithm is on the "unsmoothed" unwarped image
p_init(1, 4) = p_init(1, 4) + 0.5;
p_init(2, 4) = p_init(2, 4) + 0.5;
p_init(3, 4) = p_init(3, 4) + 0.5;

% Pick a total of n_freq_tests point offsets from pt_offsets randomly
ind         = round(size(pt_offsets, 1)*rand(n_freq_tests, 1));
ind(ind==0) =1;
pt_offsets1 = pt_offsets(ind, :);

% Scale point offsets to have required sigma
pt_offsets1 = pt_offsets1 * spatial_sigma;

% Need image to be doubles
tdata.img1 = double(tdata.img1);
tdata.img2 = double(tdata.img2);

% Space for results
results = [];

% Test counters
offset_idx = 1; 

% Convergence counters in field 1
for l = 1:length(alg_list)
    results = setfield(results, {1}, alg_list{l}, {1}, 'n_converge', 0);
end

img = tdata.img1;

% Run
while offset_idx <= n_freq_tests
    if verbose == 1
    disp(['Divergence Test: ', num2str(offset_idx)]);
    end
    % Test points: apply current point offset to target points
    test_pts = target_affine + reshape(pt_offsets1(offset_idx,:), 3, 4);
    % Solve for affine warp
    M = [template_affine; ones(1,size(template_affine,2))]' \ [test_pts; ones(1,size(test_pts,2))]';
    M = M';
    % Warp original image to get test "template" image
    tmplt = polytocuboid(tdata.img2, template_pts, M);
    % Initial error in affine points. This is not quite sqrt(mean(pt_offset(offset_idx,:) .^ 2)) due to p_init
    rms_pt_init = ComputePointError(test_pts, template_affine, p_init);
    % Run each algorithm
    
    for l=1:length(alg_list)
         if verbose == 1
            string = [alg_list{l} ' fitting...'];
            disp(string)
         end
        save matfileint.mat p_init;
        string = ['fitt = ', alg_list{l}, '(img, tmplt, p_init, n_iters, verbose, smoothing);'];
        eval(string);
        rms_pt_error = ComputePointError(test_pts, template_affine, fitt(end).warp_p);
        fitt = [];
        % Save spatial errors
        results = setfield(results, {1}, alg_list{l}, {offset_idx}, 'rms_pt_error', rms_pt_error);     
    end
    % Evaluate final spatial errors for all algorithms
    if verbose == 1
        disp('----------------------------------------------------');
        disp(['Initial spatial rms = ',num2str(rms_pt_init)]);
        disp(['Fitting results:']);
    end
    for l=1:length(alg_list)
        final_rms_pt_error = getfield(results, {1}, alg_list{l}, {offset_idx}, 'rms_pt_error');
        final_rms_pt_error = final_rms_pt_error(end);
        string = [alg_list{l}, ': final spatial rms = ', num2str(final_rms_pt_error)];
        if final_rms_pt_error < max_spatial_error
            % Update convergence counter in first field
            n_converge = getfield(results, {1}, alg_list{l}, {1}, 'n_converge');
            n_converge = n_converge + 1;
            results = setfield(results, {1}, alg_list{l}, {1}, 'n_converge', n_converge);
        else
            string = [string, '   **DIVERGED**'];
        end
        if verbose == 1
            disp(string);
        end
    end
    if verbose == 1
        disp('Display of fitting results finished, press any key to continue');
        pause; 
        disp('----------------------------------------------------');
    end
    % Always increment index into point offsets
    offset_idx = offset_idx + 1;
    if offset_idx > size(pt_offsets1, 1)+1;
        disp('Ran out of tests!');
        return;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rms_pt_error = ComputePointError(test_pts, template_affine, warp_p)
% Compute affine points rms error
% Affine for this iteration
M = build_3d_warp_a(warp_p);
% Affine points
iteration_pts = M * [template_affine; ones(1, size(template_affine, 2))];
% Error in affine points
diff_pts = test_pts - iteration_pts(1:3,:);
rms_pt_error = sqrt(mean(diff_pts(:) .^ 2));
