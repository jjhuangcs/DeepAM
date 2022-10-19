%% Find the scaling factor for IPAD thresholds
function [T,S] = threshScaleIPAD(Omega, SlT, ShT, SlV, ShV, T)
%% With the Testing and Validatin Set for Threshold Scale Selection
%%
O           = double(Omega);

%% Prediction
count = 1;
smin = 0;
step = 1e-7;
smax = step*10;

ErrorPd(count) = 1;
while(min(ErrorPd)==ErrorPd(count))
    count = count - 1;
    for s = smin:step:smax
        count = count + 1;
        AlphaL1t    = wthreshSplit(O,SlT,1e2,'s',T*s);
        Dtmp        = LQSplit(AlphaL1t,ShT,1e3,1e-10);
        AlphaL1v    = wthreshSplit(O,SlV,1e2,'s',T*s);
        Xtmp        = Dtmp*AlphaL1v;

        ErrorPd(count) = immse(Xtmp,ShV);
        PnnzPd(count) = nnz(AlphaL1v)/size(AlphaL1v,1)/size(AlphaL1v,2);
        Scale(count) = s;
    end
    
    smin = smax;
    smax = smax*10;
    step = step*10;
end

%%
[~,index]=min(ErrorPd);
S = Scale(index);
T = S*T;
end