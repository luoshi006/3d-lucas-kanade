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
verbose = 1;					% Show fitting?
smoothing = 0;                  % Apply Gaussian smoothing?
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
results = zeros(1, length(all_spc_sig), length(alg_list));

% template img
data.img1 = tdata.tmplt_img;
data.img2 = tdata.tmplt_img;
data.tmplt = tdata.tmplt;

for s=1:length(all_spc_sig)
    spatial_sigma = all_spc_sig(s);

    res = test_affine_3d(data, pt_offset, alg_list, n_iters, n_freq_tests, spatial_sigma, max_spatial_error, verbose, smoothing);
    res = cellfun(@(x) struct2cell(x), struct2cell(res), 'UniformOutput', false);
    results(1, s, :) = cell2mat(cellfun(@(x) x{1}, res, 'UniformOutput', false));
end
save(sprintf('results%d.mat', iter), 'results');