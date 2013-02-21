function create_gif(img, filename, varargin)
%% Parse input
inp = inputParser;

inp.addRequired('img', @(x)numel(size(x)) == 3);
inp.addRequired('filename', @(x)ischar(x));

inp.addOptional('DelayTime', 0.05, @(x) x >= 0 && x <= 655);
inp.addOptional('LoopCount', Inf, @(x) (x >= 0 && x < 65535) || isinf(x));

inp.parse(img, filename, varargin{:});
arg = inp.Results;
clear('inp');

%% Create gif

colormap(gray(256));

[h, w, d] = size(img);
im = zeros(h, w, 1, d);

for k=1:size(img, 3)
  im(:,:,1,k) = (img(:,:,k) / max(max(max(img(:,:,k))))) * 256;
end
imwrite(im, filename, 'DelayTime', arg.DelayTime, 'LoopCount', arg.LoopCount)