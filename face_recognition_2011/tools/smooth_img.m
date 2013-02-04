function [ simg ] = smooth_img(img)
%SMOOTH_IMG Smooth img with 5x5 Gaussian kernel

H    = fspecial('gaussian', [5 5], 2.0);
simg = imfilter(img, H, 'replicate');

end

