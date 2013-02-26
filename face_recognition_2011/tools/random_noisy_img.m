function [noisy_img, params] = random_noisy_img(tmplt, noise)

scaled = double(noise) * (max(max(max(tmplt))) - 1);
scaled(scaled == 0) = NaN;
ang = rand(1) * pi;
scaled = imrotate(scaled, ang);

noisy_img = tmplt;     

[tmplt_w, tmplt_h, tmplt_d] = size(tmplt);
[noise_w, noise_h, noise_d] = size(scaled);
w = randsample(tmplt_w - noise_w, 1);
h = randsample(tmplt_h - noise_h, 1);
d = randsample(tmplt_d - noise_d, 1);

noisy_img(w:w + noise_w - 1, h:h + noise_h - 1, d:d + noise_d - 1) = scaled;
noisy_img(isnan(noisy_img)) = tmplt(isnan(noisy_img));

params.w = w;
params.h = h;
params.d = d;
params.ang = ang;