%% Find the ratio between different thresholds
function Thre = findThreRatio(Omega, X)

sz = size(Omega,1);

A = Omega*X;
b = zeros(sz,1);
for i = 1:sz
    [~, b(i)] = fitLaplace(A(i,:));
end
clear A

Thre = 1./b;
end