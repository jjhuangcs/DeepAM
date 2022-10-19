%% Laplacian fitting
function [m, b] = fitLaplace(x)
xtmp = x;
m = 0;
b = mean(abs(xtmp-m));
end