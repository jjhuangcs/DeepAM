%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
%Images    => Cell of images used to extract patches or path to folder
%dim       => dimension of single patch
%n_patches => Total Number of patches

%% Extract training data
function [XH, XL] = trainData(Images, Patch_size, sf, stride, sigma)
XHtmp = [];
XLtmp = [];
if  ischar(Images)
    Image_list = dir(Images);
    Image_names = {};
    k = 0;
    for i=1:numel(Image_list)
        if ~Image_list(i).isdir && ~strcmp(Image_list(i).name,'Thumbs.db')
            k = k + 1;
            Image_names{k} = [Images,Image_list(i).name];
        end
    end
    for i = 1:k
        imH = imread(Image_names{i});
        imHY = imageTransformColor(imH);
        imHY = double(imHY);
        imHY = imageModcrop(imHY,sf);
        
        imLY = imageDownsample( imHY, sf, 'bicubic' );%down-sample
        
        if sigma ~=0
            imLY = double(imLY) + sigma*randn(size(imLY)); % noisy image
        end
        %%
        [Xh, Xl] = extractTrainingSetSR(imHY, imLY, Patch_size, sf, stride);
        XHtmp  = [XHtmp,Xh];
        XLtmp  = [XLtmp,Xl];
    end
end
    XH = XHtmp;
    XL = XLtmp;
end



