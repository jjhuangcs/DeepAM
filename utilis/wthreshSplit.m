%% splis the large matrix into small matrices and perform thresholding
function Y = wthreshSplit(Omega,X,splitNum,sorh,b)

[sy sx] = size(X);

for i = 1:splitNum
    tmp1 = floor(1+(i-1)*sx/splitNum);
    tmp2 = floor(i*sx/splitNum);
    tmp = wthreshVec(Omega*X(:,tmp1:tmp2),sorh,b);
    
    if i == 1
        Y = tmp;
    else
        Y = [Y,tmp];
    end
    clear tmp;
end

end