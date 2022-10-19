%% Paralle transport that transports a vector from one tangent space to the other
% (c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
% Muenchen, 2012. Contact: simon.hawe@tum.de 
function tau = parallel_transport(Phi, X, Xi, t, Norm_Xi)
if min(size(X)) == 1
    tau = Xi;
    return
end
if isempty(Phi)
    tau  = bsxfun(@times,Xi,cos(t.*Norm_Xi))-bsxfun(@times,X,(sin(t.*Norm_Xi).*Norm_Xi));
else
    tau = Phi - bsxfun(@times,...
        bsxfun(@times,X,Norm_Xi.*sin(t.*Norm_Xi)) + bsxfun(@times,Xi,(1-cos(t.*Norm_Xi))),...
        sum(Xi.*Phi)./(Norm_Xi.^2));
    tau(:,Norm_Xi==0)=Phi(:,Norm_Xi==0);
end

