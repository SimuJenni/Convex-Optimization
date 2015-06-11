
clear all;
close all;
clc;

im = im2double( rgb2gray( imread('fruits.png') ) );
% im = im2double( rgb2gray( imread('cat.png') ) );

SRfactor = 2;

% build SRfactor x SRfactor downsample operator
[M, N] = size(im);
MD = M/SRfactor;
ND = N/SRfactor;
[x,y] = meshgrid(1:N,1:M);
cols = (x(:)-1)*SRfactor*MD+y(:);
rows = (floor((x(:)-1)/SRfactor))*MD+floor((y(:)-1)/SRfactor)+1;
vals = 1/SRfactor^2*ones(M*N,1);
D = sparse(rows,cols,vals,MD*ND,M*N);

% create input image
g = reshape(D*im(:),M/SRfactor,N/SRfactor);

lambda = 90;

uG = superresolution_SIMONJENNI(g,D,lambda);

figure;
disp = [uG, 0*im; ...
        im, imresize(g,[M N],'nearest') ];
imshow(disp);

sum(sum((uG-im).^2))
