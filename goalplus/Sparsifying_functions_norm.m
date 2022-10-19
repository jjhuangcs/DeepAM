%Sparsifying_functions Computes various sparsity promoting functions.
%   Type. Type defines the type of sparsifying function to be used. Possible
%   values are:
%   'LogAbs' Sum_i log(1+mu*|X|) or its derviative depending on Purpose
%   'LogSquare' Sum_i log(1+mu*X^2) or its derviative depending on Purpose
%   'PNormAbs' Sum_i (|X|+mu)^p or its derviative depending on Purpose
%   'PNormSquare' Sum_i (|X|^2+mu)^p/2 or its derviative depending on Purpose
%   'AtanAbs' Sum_i atan(mu*|X|) or its derviative depending on Purpose
%   'AtanSquare' Sum_i atan(mu*|X|.^2) or its derviative depending on Purpose7
%   Purpose
%   'Derivative' the derivative of the above mentioned functions is
%   computed
%   'Otherwise' the function value is computed
%   X. X is the Matrix/Vector whose sparsity should be computed
%   Exponents. If the variable Exponents is not empty, it is assumed that
%   the mixed q-Type sparsity measure has to be computed, and q =
%   Exponents(1). For the p-norms p = Exponents(2)
%   n. If n is not empty, the vectors in X are assumed to be in block form
%   and each block contains n elements. This simplyfies the computation of
%   the deriviative for mixed q-Type norms.
%   q_weight. Is for computing the derivative of mixed q-Type norms, if the
%   function value of the corresponding function has been computed
%   previously


function [f,q_weight] = Sparsifying_functions_norm(Type, Purpose, X, Exponents, mu, n, q_weight)
% Depending on n, X is assorted in blocks or the signals are stacked
if nargin == 5 || isempty(n)
    sz = size(X);
    n = prod(sz(1:end-1));
    %     n = size(X,1);
    %     if size(X,3) ~= 1
    %         n = n*size(X,2);
    %     end
end

if ~isempty(Exponents)
    q = Exponents(1);
else
    q = 0;
end

if nargin == 6
    q_weight = 1;
end

% Check whether to compute derivatives or function values
derivative = strcmpi(Purpose,'Derivative');
Comp_qweight = (q ~= 0  & ~(nargin == 7));
switch Type
    case 'LogAbs'
        % log(1+mu*|x|)
        if ~derivative || Comp_qweight
            q_weight = sum(reshape(log(1+mu*abs(X)),n,[]));
        end
        if derivative
            f = mu*sign(X)./(1+mu*abs(X));
        end
    case 'LogSquare'
        % log(1+mu*x^2)
        if ~derivative || Comp_qweight
            q_weight = sum(reshape(log(1+mu*X.^2),n,[]))/log(1+mu);
        end
        if derivative || Comp_qweight
            f = mu*2*X./(1+mu*X.^2)/log(1+mu);
        end
        
    case 'LogSquare2'
        % log(1+mu*x^2)
        if ~derivative || Comp_qweight
            q_weight = sum(reshape(log(1+mu*sqrt(X.^2+1e-6)),n,[]));
        end
        if derivative || Comp_qweight
            ssq = sqrt(X.^2+1e-6);
            f = mu*X./((1+mu*ssq).*ssq);
        end
        
    case 'ExpSquare'
        % 1 - exp(-x^2/(2mu))
        if ~derivative || Comp_qweight
            %q_weight = sum(reshape(1-exp(-X.^2/(2*mu)),n,[]));
            w = 1/(sqrt(pi)*mu);
            q_weight = sum(reshape(w-w*exp(-X.^2/(mu^2)),n,[]));
        end
        if derivative || Comp_qweight
            w = 1/(sqrt(pi)*mu);
            %f = X./mu.*exp(-X.^2/(2*mu));
            f = 2*w*X./mu.*exp(-X.^2/(mu^2));
        end
        
    case 'ExpAbs'
        % 1 - exp(-abs(x)/(2mu))
        if ~derivative || Comp_qweight
            q_weight = sum(reshape(1-exp(-abs(X)/(mu)),n,[]));
        end
        if derivative || Comp_qweight
            f = sign(X)./mu.*exp(-abs(X)/(mu));
        end
    case 'PNormAbs'
        p = Exponents(2);
        % (|x|+mu)^p
        if ~derivative
            if p ~= 1
                q_weight = sum(reshape(exp(p*log(abs(X)+mu)),n,[]));
            else
                q_weight = sum(reshape(abs(X),n,[]));%l1 norm
            end
        end
        if derivative || Comp_qweight
            if p~=1
                f = p*sign(X).*exp((p-1)*log(abs(X)+mu));
            else
                f = sign(X);
            end
        end
    case 'PNormSquare'
        p = Exponents(2);
        % (|x|^2+mu)^p/2
        if ~derivative
            q_weight = sum(reshape(exp(p/2*log(X.^2+mu)),n,[]));
        end
        if derivative || Comp_qweight
            f = p/2*X.*exp((p/2-1)*log(X.^2+mu));
        end
    case 'AtanAbs'
        % atan(mu*|x|)
        if ~derivative || Comp_qweight
            q_weight = sum(reshape(atan(abs(X)*mu),n,[]));
        end
        if derivative
            f = mu*sign(X)./(1+(mu*abs(X)).^2);
        end
    case 'AtanSquare'
        % atan(mu*|x|.^2)
        if ~derivative || Comp_qweight
            q_weight = sum(reshape(atan(X.^2*mu),n,[]));
        end
        if derivative
            f = 2*mu*X./(1+(mu*X).^4);
        end
end

if derivative
    f = reshape(bsxfun(@times,reshape(f,n,[]),q_weight.^(q-1)),size(X));
end

if ~derivative
    if q ~=0
        %f = 1/(2*q)*sum((q_weight.^q));
        f = 1/(q)*sum((q_weight.^q));
        %q_weight = q_weight.^(q-1);
    else
        f = sum(q_weight);
        q_weight = 1;
    end
end