%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
%Images    => Cell of images used to extract patches or path to folder
%dim       => dimension of single patch
%n_patches => Total Number of patches

%% Modified by Junjie 2017-10-09
function [XH, XL, imLY] = completeTestingSetExtractSR(Images, Patch_size, sf, stride, sigma)
XHtmp = [];
XLtmp = [];
if  ischar(Images)
    imH = imread(Images);
    imHY = imageTransformColor(imH);
    imHY = double(imHY);
    imHY = imageModcrop(imHY,sf);
    imLY = imageDownsample( imHY, sf, 'bicubic' );%down-sample
    %%
    if sigma ~= 0
        imLY = double(imLY) + sigma*randn(size(imLY)); %noisy image
    end
    %%
    [Xh, Xl] = extractTrainingSetSR(imHY, imLY, Patch_size,sf, stride);
    XHtmp  = [XHtmp,Xh];
    XLtmp  = [XLtmp,Xl];
end
    XH = XHtmp;
    XL = XLtmp;
end



