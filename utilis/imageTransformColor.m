function [imY, imCB, imCR] = imageTransformColor( im )
% imageTransformColor Transform an image into the YCBCR color space.
% 
% Transforms a given image im into the YCBCR color space and returns the
% separate components (Y, CB, and CR). In case the given image im is a
% gray-scale image, the CB and CR components are empty. 
% 
% INPUTS
%  im           - [imh x imw x nchs] an image
% 
% OUTPUTS
%  imY          - [imh x imw] the Y channel of the image
%  imCB         - [imh x imw] the CB channel of the image
%  imCR         - [imh x imw] the CR channel of the image
% 
% Code adapted from [1]
% 
% References:
% [1] R. Timofte, V. De Smet, L. van Gool. Anchored Neighborhood Regression 
% for Fast Example-Based Super- Resolution. ICCV 2013. 

if size(im,3) == 3
  im = rgb2ycbcr(im); imY=im2single(im(:,:,1));
  imCB = im2single(im(:,:,2)); imCR = im2single(im(:,:,3));
else
  imY = im2single(im); imCB = []; imCR = [];
end
end