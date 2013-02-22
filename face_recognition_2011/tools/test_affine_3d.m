function results = test_affine_3d(tdata, pt_offsets, alg_list, n_iters, n_tests, n_freq_tests, spatial_sigma, image_pixel_sigma, tmplt_pixel_sigma, max_spatial_error, verbose, smoothing)
% TEST_AFFINE - Test affine algorithms
%
% See also: run_affine
%
% tdata has two fields:
%   tdata.img       unit8 greyscale image
%   tdata.tmplt     [x_start y_start x_end y_end] template rectangle corners
%                   matlab coordinates (1, 1) is top-left, start < end.
%
% pt_offsets is a N x 6 matrix of (x1, y1, x2, y2, x3, y3) deltas for each
% of the three affine control points. Zero mean, unit variance, 
% e.g. pt_offsets = randn(N, 6);
%
% alg_list is a cellstr list of algorithms to run, e.g.:
% alg_list = {'affine_fa' 'affine_fc' 'affine_ia' 'affine_ic'};
%
% The "template" is created by cutting out a distorted version of tdata.tmplt.
% The "target" is the image defined by tdata.tmplt.

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: test_affine.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

if nargin<11 error('Not enough input arguments'); end

minX = tdata.tmplt(1);
minY = tdata.tmplt(2);
minZ = tdata.tmplt(3);
maxX = tdata.tmplt(4);
maxY = tdata.tmplt(5);
maxZ = tdata.tmplt(6);

% Target corner points (i.e. correct answer)
target_pts = [minX minY minZ; % Bottom rectangle
              minX maxY minZ;
              maxX maxY minZ;
              maxX minY minZ;
              minX minY maxZ; % Top rectangle
              minX maxY maxZ;
              maxX maxY maxZ;
              maxX minY maxZ]';

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

% Scale point offsets to have required sigma
pt_offsets = pt_offsets * spatial_sigma;

% Need image to be doubles
tdata.img = double(tdata.img);

% Space for results
results = [];

% Test counters
go = 1; offset_idx = 1; all_alg_converged = 0;

% Convergence counters in field 1
for l=1:length(alg_list)
	results = setfield(results, {1}, alg_list{l}, {1}, 'n_diverge', 0);
	results = setfield(results, {1}, alg_list{l}, {1}, 'n_converge', 0);
end

% Index of successfull convergence tests
results.all_converged_idx = [];

% Might not be doing any convergence tests at all
if n_tests > 0
	cnvg_testing = 1;
else
	cnvg_testing = 0;
	cnvg_result = [];
end

% Run
while go
	if cnvg_testing
		disp(['Convergence Test: ', num2str(all_alg_converged + 1), ' (of total ', num2str(offset_idx), ')']);
	else
		disp(['Divergence Test: ', num2str(offset_idx)]);
	end 
	
	% Test points: apply current point offset to target points
	test_pts = target_affine + reshape(pt_offsets(offset_idx,:), 3, 4);
		
	% Solve for affine warp
	M = [template_affine; ones(1,size(template_affine,2))]' \ [test_pts; ones(1,size(test_pts,2))]';
	M = M';
	
	% Warp original image to get test "template" image
	tmplt_img = polytocuboid(tdata.img, template_pts, M);
	
	% Add noise to template
	if tmplt_pixel_sigma > 0
		tmplt_img = tmplt_img + (randn(size(tmplt_img)) * tmplt_pixel_sigma);
	end
	
	% Add noise to image
	if image_pixel_sigma > 0
		noisy_img = tdata.img + (randn(size(tdata.img)) * image_pixel_sigma);
	else
		noisy_img = tdata.img;
	end
	
	% Initial error in affine points. This is not quite sqrt(mean(pt_offset(offset_idx,:) .^ 2)) due to p_init
	rms_pt_init = ComputePointError(test_pts, template_affine, p_init);
	
	% Run each algorithm
	for l=1:length(alg_list)
        fitp = [];
		string = ['tic; fitp = ', alg_list{l}, '(noisy_img, tmplt_img, p_init, n_iters, verbose, smoothing); t = toc;'];
		eval(string);
		
		% Evaluate point spatial error for each iteration for convergence tests
		if cnvg_testing
			rms_pt_error = zeros(length(fitp), 1);
			for f=1:length(fitp)
				rms_pt_error(f) = ComputePointError(test_pts, template_affine, fitp(f).warp_p);
			end
					
		% Only need final spatial rms for divergence tests. Don't need fitting results
		else
			rms_pt_error = ComputePointError(test_pts, template_affine, fitp(end).warp_p);
		end	
		
		% Save spatial errors
		results = setfield(results, {1}, alg_list{l}, {offset_idx}, 'rms_pt_error', rms_pt_error);
			
		% Save full fitting results
		results = setfield(results, {1}, alg_list{l}, {offset_idx}, 'fit', fitp);
			
		% Save fitting time
		results = setfield(results, {1}, alg_list{l}, {offset_idx}, 'time', t);
	end
	
	% Evaluate final spatial errors for all algorithms
	disp('----------------------------------------------------');
	disp(['Initial spatial rms = ',num2str(rms_pt_init)]);
	
	all_cnvg_check = 1;
	
	for l=1:length(alg_list)
		final_rms_pt_error = getfield(results, {1}, alg_list{l}, {offset_idx}, 'rms_pt_error');
		final_rms_pt_error = final_rms_pt_error(end);

		string = [alg_list{l}, ': final spatial rms = ', num2str(final_rms_pt_error)];
		
		% Convergence test
		if final_rms_pt_error > max_spatial_error
			string = [string, '   **DIVERGED**'];
			
			% Update divergence counter in first field
			n_diverge = getfield(results, {1}, alg_list{l}, {1}, 'n_diverge');
			n_diverge = n_diverge + 1;
			results = setfield(results, {1}, alg_list{l}, {1}, 'n_diverge', n_diverge);
			
			% One or more algorithms diverged
			all_cnvg_check = 0;
			
		else
			% Update convergence counter in first field
			n_converge = getfield(results, {1}, alg_list{l}, {1}, 'n_converge');
			n_converge = n_converge + 1;
			results = setfield(results, {1}, alg_list{l}, {1}, 'n_converge', n_converge);
        end
        
%        % Display error
%         clf;
%         error_img = abs(tmplt_img - warp_3d_a(noisy_img, fitp(end).warp_p, template_pts));
%         error_img(error_img == 0) = NaN;
%         PATCH_3Darray(error_img, 'col', 'barS');
%         view([-12 8]);
		
		disp(string);
	end
	disp('----------------------------------------------------');
	
	% Everything converged?
	if all_cnvg_check
		all_alg_converged = all_alg_converged + 1;
		
		% Update index of fully converged test results
		results.all_converged_idx = [results.all_converged_idx; offset_idx];
	end
	
	% Finished convergence testing?
	if all_alg_converged >= n_tests
		cnvg_testing = 0;
		
		% Extra tests for divergence?
		if n_freq_tests > n_tests
			if offset_idx >= n_freq_tests
				go = 0;
			end
			
		% Or just stop
		else
			go = 0;
		end
	end
	
	% Always increment index into point offsets
	offset_idx = offset_idx + 1;
	if offset_idx > size(pt_offsets, 1)
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
