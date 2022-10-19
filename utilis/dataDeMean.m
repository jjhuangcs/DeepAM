%% remove the mean of the patches
function [out] = dataDeMean(in)

out = deMeanOperator(sqrt(size(in,1)))*in;    

end