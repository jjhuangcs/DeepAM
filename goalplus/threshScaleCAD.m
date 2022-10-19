%% Find the scaling factor for the CAD thresholds
function [T,S] = threshScaleCAD(Omega, SlT, ShT, SlV, ShV, Tg, Tl)
% With the Testing and Validatin Set for Threshold Scale Selection

% Fix the threhsold for the local atoms, while finding threshold scale for
% the global atoms.
%%
O               = double(Omega);
ratio           = (size(Tg,1) - size(Tl,1))/size(Tg,1);

%% Prediction
ifBreak         = 0;
count           = 1;
smin            = 0;
step            = 1e-2;
smax            = step*10;

ErrorPd(count)  = 1;
PnnzPd(count)   = 1;

while((PnnzPd(count)>2*ratio)||(ErrorPd(count)<=ErrorPd(1)))
    count = count - 1;
    for s = smin:step:smax
        count = count + 1;
        Ttmp = Tg;
        Ttmp(end - size(Tl,1) + 1:end) = Tl*s;
        AlphaL1t    = wthreshSplit(O,SlT,1e2,'s',Ttmp);
        Dtmp        = LQSplit(AlphaL1t,ShT,1000,1e-10);
        
        AlphaL1v    = wthreshSplit(O,SlV,1e2,'s',Ttmp);
        Xtmp        = Dtmp*AlphaL1v;
        
        ErrorPd(count)  = immse(Xtmp,ShV);
        PnnzPd(count)   = nnz(AlphaL1v)/size(AlphaL1v,1)/size(AlphaL1v,2);
        Scale(count)    = s;
        
        if (PnnzPd(count)==ratio)
            ifBreak = 1;
            break;
        end
    end
    if ifBreak==1 || PnnzPd(count)==PnnzPd(count-1)
        break;
    end
    smin = smax;
    smax = smax*10;
    step = step*10;
end

%%
[minErr,index]=min(ErrorPd);
S = Scale(index);
T = Tg;
T(end - size(Tl,1) + 1:end) = Tl*S;

S=minErr;
end