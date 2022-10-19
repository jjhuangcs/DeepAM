% (c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
% Muenchen, 2012. Contact: simon.hawe@tum.de
%% Exponentiall mapping 
function D = exp_mapping(D, dX, t, Norm_dx,sel)

if min(size(D))==1
    D = D+t*dX;
end

D(:,sel) = bsxfun(@times,D(:,sel),cos(t.*Norm_dx(:,sel)))+...
               bsxfun(@times,dX(:,sel),(sin(t.*Norm_dx(:,sel))./Norm_dx(:,sel)));

D =  bsxfun(@times,D, 1./sqrt(sum(D.^2))); %Normalization due to numerical things
					


