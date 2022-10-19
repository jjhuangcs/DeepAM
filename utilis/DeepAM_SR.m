function [imSR,imBic,psnrSR]=DeepAM_SR(p_sz0, sf, Images, Omega, T, sd, sigma)
ifrm        = 1;% remove mean
ifSubMat    = 1;% sub-matrix
%%
[XHt, XLt, imBic] =  completeTestingSetExtractSR(Images, p_sz0, sf, sd, sigma);

if ifrm ==1
    Slt     = dataDeMean(XLt);
    SlMean  = mean(XLt); 
else
    Slt = XLt;
end

[~,szc] = size(T);
for i=1:szc
    Slt = wthreshSplit(Omega{i},Slt,1,'s',T{i});
end
Y   = Omega{szc+1}*Slt;
if ifrm==1
    Y   = Y + ones(size(Y,1),1)*SlMean;% Add the mean back
end
if ifSubMat == 1
    pszSub  = 8;
    subMat  = extractCenterPatch(Y, pszSub);
    Y       = subMat*Y;
end
%% Image SR
imH     = imread(Images);
imH     = imageModcrop(imH, sf);
img_size= size(imH);

gridH   = getSamplingGrid(size(imH),sf*[p_sz0,p_sz0],sf*[p_sz0,p_sz0]-sd*[1,1]*sf,[1,1]*sf,1);

if ifSubMat == 1
    imSR    = uint8(255*overlap_add2( Y, img_size, gridH ));
else
    imSR    = uint8(255*overlap_add( Y, img_size, gridH ));
end

% shave
imH = shave(imH,p_sz0*[sf,sf]);   
imSR = shave(imSR,p_sz0*[sf,sf]); 

% psnr
psnrSR  = csnr_index(double(imH),double((imSR)),0,0);

end