function [noisy_img] = random_img_per_slice(tdata, w, h)

tmplt = tdata.tmplt_img;

box = tdata.tmplt;
min_w = box(1);
min_h = box(2);
min_d = box(3);
max_w = box(4);
max_h = box(5);
max_d = box(6);

noisy_img = tmplt;

for i = min_d:max_d - 1
    minx = randsample(1:(max_w - w), 1);
    miny = randsample(1:(max_h - h), 1);
    z = randsample(size(tmplt, 3), 1);

    slice = tmplt(miny:miny + h - 1, minx:minx + w - 1, z);
    
    rand_w = randsample(min_w:(max_w - w), 1);
    rand_h = randsample(min_h:(max_h - h), 1);
    
    noisy_img(rand_h:rand_h + h - 1, rand_w:rand_w + w - 1, i) = slice;
end