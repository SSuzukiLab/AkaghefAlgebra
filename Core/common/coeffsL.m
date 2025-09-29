function [C,T]=coeffsL(p,var,varargin)
% coeffs for Laurent polynomials
    if nargin==1
        var = symvar(p); % if not provided, use variables from p
    end
    d=zeros(size(var));
    for i=1:length(var)
        tmp=series(p,var(i),0,Order=1);
        % degree of polynomial
        if tmp==0,continue; end
        try
            d(i)=-polynomialDegree(tmp,var(i));
        catch
            d(i)=polynomialDegree(1/tmp,var(i));
        end
    end
    
    p2= simplify(p * prod(var.^d));
    [C,T]=coeffs(p2,var);
    T= T ./ prod(var.^d);
    if nargin==3&&varargin{1}=="int"
        assert(isscalar(var))
        assume(var>0)
        T=double(simplify(log(T)/log(var)));
    end
end