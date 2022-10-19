%% Perform least squares in a batch manner
function D = LQSplit(X,Y,splitNum,lambda)

sx = size(X,2);

covaMat = zeros(size(Y,1),size(X,1));
autoMat = zeros(size(X,1),size(X,1));


for i = 1:splitNum
    tmp1 = floor(1+(i-1)*sx/splitNum);
    tmp2 = floor(i*sx/splitNum);
    
    TmpX = X(:,tmp1:tmp2);
    TmpY = Y(:,tmp1:tmp2);
    
    autoMat = autoMat + TmpX*TmpX';
    covaMat = covaMat + TmpY*TmpX';
end

D = covaMat*pinv(autoMat + lambda*eye(size(autoMat)));

end