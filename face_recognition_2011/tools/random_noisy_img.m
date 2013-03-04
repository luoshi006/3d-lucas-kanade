function [noisy_img, params] = random_noisy_img(tdata, s)

tmplt = tdata.tmplt_img;

box = tdata.tmplt;
min_w = box(1);
min_h = box(2);
min_d = box(3);
max_w = box(4);
max_h = box(5);
max_d = box(6);

tmplt_width = max_w - min_w - 1;
tmplt_height = ((max_h - min_h - 1) / 2) * s;
tmplt_depth = max_d - min_d - 1;

minx = randsample(1:(max_w - tmplt_width), 1);
miny = randsample(1:(max_h - tmplt_height), 1);
minz = randsample(1:(max_d - tmplt_depth), 1);

noise = tmplt(miny:(miny + tmplt_height), minx:(minx + tmplt_width), minz:(minz + tmplt_depth));
noise(noise == 0) = NaN;

noisy_img = tmplt;

[noise_h, noise_w, noise_d] = size(noise);
% w = randsample(min_w:(max_w - noise_w), 1);
h = randsample(min_h:(max_h - noise_h), 1);
% d = randsample(min_d:(max_d - noise_d), 1);

noisy_img(h:h + noise_h - 1, min_w:min_w + noise_w - 1, min_d:min_d + noise_d - 1) = noise;
noisy_img(isnan(noisy_img)) = tmplt(isnan(noisy_img));

params.w = min_w;
params.h = h;
params.d = min_d;