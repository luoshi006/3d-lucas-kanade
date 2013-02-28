function [noisy_img, params] = random_noisy_img(tdata, scale)

tmplt = tdata.tmplt_img;
noise = tdata.noise_img;
box = tdata.tmplt;

scaled = double(noise) * (max(max(max(tmplt))) - 1);
scaled = imresize(scaled, scale, 'bilinear');
ang = rand(1) * pi;
scaled = imrotate(scaled, ang);
scaled(scaled == 0) = NaN;

noisy_img = tmplt;     

min_w = box(1);
min_h = box(2);
min_d = box(3);
tmplt_w = box(4) - min_w;
tmplt_h = box(5) - min_h;
tmplt_d = box(6) - min_d;
[noise_w, noise_h, noise_d] = size(scaled);
w = randsample(min_w:tmplt_w - noise_w, 1);
h = randsample(min_h:tmplt_h - noise_h, 1);
d = randsample(min_d:tmplt_d - noise_d, 1);

noisy_img(w:w + noise_w - 1, h:h + noise_h - 1, d:d + noise_d - 1) = scaled;
noisy_img(isnan(noisy_img)) = tmplt(isnan(noisy_img));

params.w = w;
params.h = h;
params.d = d;
params.ang = ang;