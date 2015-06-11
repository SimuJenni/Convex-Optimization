clear all;
close all;
clc;

im = im2single( rgb2gray( imread('fruits.png') ) );

% noisy image
sigma = 0.1;
noise = sigma * randn(size(im));
g = max( min( noise + im, 1 ), 0 );

lambda = 10;

uG = denoising_Dummy(g,lambda);

disp = [ uG, 0*im; ...
         im, g ];

figure;
imshow(disp);
