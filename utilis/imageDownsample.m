function imL = imageDownsample( imH, sf, kernel )
% imageDownsample Down-scales a given image.
%
% Given a high-resolution image imH, this function downscales it to a
% low-resolution image imL by the factor sf. The factor defines the amount
% of DOWN-sampling, i.e., sf=2 means that the image size is halved. The
% options of the function include different down-sampling kernels and the
% corresponding settings. Currently, only the bicubic kernel is
% implemented.

switch kernel
    case 'bicubic'
        imL = imresize(imH,1/sf,'bicubic');
    case 'Gaussian'
        error('Not implemented yet!');
    case 'nearest'
        imL = imresize(imH,1/sf,'nearest');
    case 'bilinear'
        imL = imresize(imH,1/sf,'bilinear');   
    otherwise
        imL = imresize(imH,1/sf,kernel);
end
end