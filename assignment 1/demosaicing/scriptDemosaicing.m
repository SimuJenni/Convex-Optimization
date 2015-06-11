
clear all;
close all;
clc;

im = im2single(imread('fruits.png'));

% masks
[m,n,~] = size(im);
blue = zeros(m,n);
blue(1:2:end,1:2:end) = 1;
red = zeros(m,n);
red(2:2:end,2:2:end) = 1;
green = ones(m,n);
green(1:2:end,1:2:end) = 0;
green(2:2:end,2:2:end) = 0;

% mosaic image
g = zeros(size(im));
g(:,:,1) = red .* im(:,:,1);
g(:,:,2) = green .* im(:,:,2);
g(:,:,3) = blue .* im(:,:,3);
gdisp = g;
g = sum(g,3);

lambda = 1000;

uG = demosaicing_Dummy(g,lambda);

figure;
disp = [uG, 0*im; ...
        im, gdisp];
imshow(disp);

