%% Initialize the dictionary
function [Omega] = initDict(Operator_type,S,PW,lift)
% PW        PatchWidth

Omega = [];
for i=1:numel(Operator_type)
    % Selecting some type of initial analysis operator
    switch Operator_type{i}
        case 'TV'
            O = create_tvmat(PW);
        case 'SVD'
            [O,~]=eig(S*S');
            O = create_lifting(lift,O);
        case 'SVDN'
            [O,vals]=eig(S*S');
            O = bsxfun(@times,O,1./diag(vals));
            O = bsxfun(@times,O,1./sqrt(sum(O.^2,1)));
            O = create_lifting(lift,O);
        case 'DCT'
            O = kron(dctmtx(PW),dctmtx(PW));
            O = create_lifting(lift,O);
        case 'RAND'
            O = randn(round(lift*size(S,1)),size(S,1));

        case 'ELAD'
            O = initOmega(S,round(lift*size(S,1)));
        case 'LOAD'
            O = load(sprintf('%dx%dSquared.mat',PW,PW));
            O = create_lifting(lift,O.Omega);
        case 'NONE'
            O = Learn_para.Omega;
    end
    Omega = [Omega;O];
end

Omega       = bsxfun(@times,Omega,1./sqrt(sum(Omega.^2,2)));