%WTHRESH Perform soft or hard thresholding. 
%   Y = WTHRESH(X,SORH,T) returns soft (if SORH = 's')
%   or hard (if SORH = 'h') T-thresholding  of the input 
%   vector or matrix X. T is the threshold value.
%
%   Y = WTHRESH(X,'s',T) returns Y = SIGN(X).(|X|-T)+, soft 
%   thresholding is shrinkage.
%
%   Y = WTHRESH(X,'h',T) returns Y = X.1_(|X|>T), hard
%   thresholding is cruder.
%
%   See also WDEN, WDENCMP, WPDENCMP.
 
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 24-Jul-2007.
%   Copyright 1995-2007 The MathWorks, Inc.
function y = wthreshVec(x,sorh,b)
 
switch sorh
  case 's'
    tmp = (abs(x)-b*ones(1,size(x,2)));
    tmp = (tmp+abs(tmp))/2;
    y   = sign(x).*tmp ;
 
  case 'h'
    tmp = (abs(x)>b*ones(1,size(x,2)));
    y   = x.*tmp;
  otherwise
    error(message('Wavelet:FunctionArgVal:Invalid_ArgVal'))
end

