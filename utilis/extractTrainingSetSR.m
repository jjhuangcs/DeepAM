%% This function extracts n training samples of size sz from Signal S.
% The input S can be either a 1-D signal or 2-D signal (image). For a 1-D
% signal n training sequences having sz samples will be extracted from
% random positions. For a 2-D Signal n training patches of size
% [sz(1),sz(2)] or [sz,sz] will be extracted from n random positions.
% (c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
% Muenchen, 2012. Contact: simon.hawe@tum.de

%% Modified by Junjie 2017-11-02
function [Xh, Xl] = extractTrainingSetSR(imH,imL, sz,sf,stride)

if isscalar(sz) && min(size(imH)) > 1
    szH = [sz,sz]*sf;
    szL = szH/sf;
end

% Compute one grid for all filters
gridH = getSamplingGrid(size(imH),szH,szH-stride*[1,1]*sf,[1,1]*sf,1);
fH = imH(gridH); Xh=reshape(fH,[size(fH,1)*size(fH,2),size(fH,3)]);

gridL = getSamplingGrid(size(imL),szL,szL-stride*[1,1],[1,1],1);
fL = imL(gridL); Xl=reshape(fL,[size(fL,1)*size(fL,2),size(fL,3)]);



