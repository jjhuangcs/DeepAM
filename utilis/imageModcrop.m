function imCrop = imageModcrop( im, modfactor )
% imageModcrop Crops an image to be divideable without remainder.
%
% Given an image im and a modulo-factor modfactor, this function crops the
% image such that it is exactly dividable by modfactor, i.e., without any
% remainder. 
% 
% INPUTS
%  im           - [imh x imw x nchns] an image
%  modfactor    - the modulo factor
% 
% OUTPUTS
%  imCrop       - the cropped image
% 
% Code adapted from [1]
% 
% References:
% [1] R. Timofte, V. De Smet, L. van Gool. Anchored Neighborhood Regression 
% for Fast Example-Based Super- Resolution. ICCV 2013. 

[imh,imw,~]=size(im); imh=imh-mod(imh,modfactor); imw=imw-mod(imw,modfactor);
imCrop = im(1:imh,1:imw,:);
end