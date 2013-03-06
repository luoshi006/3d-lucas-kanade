function [ simg ] = smooth_img(img, varargin)
%% Parse input
inp = inputParser;

inp.addRequired('img', @(x) numel(size(x)) == 3);
inp.addOptional('siz', [5 5 5], @(x) x > 0);
inp.addOptional('sigma', 2.0, @(x) x > 0);

inp.parse(img, varargin{:});
arg = inp.Results;

if length(arg.siz) == 1
    arg.siz = repmat(arg.siz, 1, 3);
end

clear('inp');

%%

x = -ceil((arg.siz - 1) / 2):ceil((arg.siz - 1) / 2);
H = exp(-(x.^2 / (2 * arg.sigma ^ 2)));
H = H / sum(H(:));
    
Hx = reshape(H, [length(H) 1 1]);
Hy = reshape(H, [1 length(H) 1]);
Hz = reshape(H, [1 1 length(H)]);
simg = imfilter(imfilter(imfilter(img, Hx, 'replicate'), Hy, 'replicate'), Hz, 'replicate');
end

