function results = run_affine_3d(data_name)
% Run an affine perturbation test
%
% e.g. run_affine_3d('takeo');
%
% You should edit this file to define experiment parameters!

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: run_affine.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

if nargin<1 error('Not enough input arguments'); end

% List of algorithms to run
alg_list = get_all_files('methods', 'affine(_[\w]+)?_ic([_A-Za-z]+)?_3d\.(p|m)');
alg_list = cellfun(@(x) x(1:length(x)-2), alg_list, 'UniformOutput', false);

% Test parameters
verbose = 0;					% Show fitting?
smoothing = 0;                  % Apply Gaussian smoothing?
scale = 0.1;                    % % increase in size increase every iteration
n_iters = 30;					% Number of gradient descent iterations
n_freq_tests = 100;				% Number of frequency of convergence tests
max_spatial_error = 1.4;        % Max location error for deciding convergence

all_spc_sig = (1:10);		    % All spatial sigmas

% Should not need to modify anything below --------------------------------

% tdata - the image and initial template
tdata = load(['data/', data_name]);

% pt_offset - precomputed random point offsets
pt_offset = load('data/affine_pt_offset_3d');
pt_offset = pt_offset.pt_offset;
num_of_scales = 10;
results = zeros(num_of_scales, length(all_spc_sig), length(alg_list));

% Run tests
for iter = 1:num_of_scales
    % template img
    data.img2 = tdata.tmplt_img;
    data.tmplt = tdata.tmplt;

    s = 1 + ((iter - 1) * scale);
    [data.img1, params] = random_noisy_img(tdata, s);
    create_gif(data.img1, sprintf('%dx%dx%d - (%f,%f,%f) (%.2f).gif', params.w,  params.h,  params.d,  params.angx, params.angy, params.angz, s));

    % Matrix S for Gabor-Fourier method, thanx to Peter Kovesi's Gabor Filters, http://www.csse.uwa.edu.au/~pk/
%         temp = ones(tdata.tmplt(4)-tdata.tmplt(2)+1, tdata.tmplt(3)-tdata.tmplt(1)+1);
%         num_of_scales = 32; num_of_or = 32;
%         [EO, BP, S] = gaborconvolve(temp, num_of_scales, num_of_or , 3, 2, 0.65);
%         save S.mat S;

    % Run tests
    parfor s=1:length(all_spc_sig)
        spatial_sigma = all_spc_sig(s);

        fprintf('Scale: %d - Sigma: %.2f \n', iter, spatial_sigma);

        res = test_affine_3d(data, pt_offset, alg_list, n_iters, n_freq_tests, spatial_sigma, max_spatial_error, verbose, smoothing);
        res = cellfun(@(x) struct2cell(x), struct2cell(res), 'UniformOutput', false);
        results(iter, s, :) = cell2mat(cellfun(@(x) x{1}, res, 'UniformOutput', false));
    end
    save(sprintf('results%d.mat', iter), 'results');
end