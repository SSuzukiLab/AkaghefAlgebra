classdef qNumS
    methods(Static)
        function out=n(q,x)
            out=(q^x-q^(-x))/(q-q^(-1));
        end
        function out=tau(q)
            out=q-q^-1;
        end
        function ret=pch(q,x,n)
            if n==0
                ret=1;
            else
                ret=(1-x*q^(n-1))*qNumS.pch(x,n-1,q);
            end
        end
        function out=fac(q,n)
            out=sym(1);
            for kk=2:n
                out=out*qNumS.n(q,kk);
            end
        end
        function out=exp(q,X,n,inv)
            arguments
                q
                X
                n
                inv (1,1) logical=false
            end
            if inv
                q=q^-1;
            end
            if nargin==2
                n=5;
            end
            out=0;
            for k=0:n
                out=out+(q^(k*(k-1)/2)/qNumS.fac(q,k))*X^k;
            end
        end
        function ret=binomial(q,n,m)
            if ~isa(n,'sym')&&n<m
                ret=0;
                return
            else
                ret=1;
            end
            for ii=0:m-1
                ret=qNumS.n(q,n-ii)/qNumS.n(q,m-ii)*ret;
            end
        end
        function ret=factor(q,arg)
            F=factor(sym(arg));
            N=numel(F);
            T=table(nan(N,1),F.',VariableNames={'type','value'});
            ret=T;
            % for i=1:height(T)
            %     x=T{i,'value'};
            %     if isAlways(isfinite(x))
            %         t=1;
            %     elseif isAlways(x==q)||isAlways(x==1/q)
            %         t=2;
            %     else
            %         [c,t] = coeffs(___)
            %     end
            % end
        end
        function ret=test(arg)
            ret=sym('q');
        end
    end
end