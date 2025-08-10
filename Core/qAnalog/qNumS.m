classdef qNumS
    properties
        q=sym('q')
    end
    methods
        function qN=qNumS(q)
            if nargin==0, q=sym('q'); end
            qN.q=q;
        end
        function out=n(qN,x)
            q=qN.q;
            out=(q^x-q^(-x))/(q-q^(-1));
        end
        function out=T(qN)
            q=qN.q;
            out=q-q^-1;
        end
        function ret=pch(qN,x,n,Q)
            if n==0
                ret=1;
            else
                if nargin==3, Q=qN.q; end
                ret=(1-x*Q^(n-1))*qN.pch(x,n-1,Q);
            end
        end
        % function out=q(obj)
        %     syms q
        %     out=q;
        % end
        function out=fac(obj,n)
            out=sym(1);
            for kk=2:n
                out=out*obj.n(kk);
            end
            % out=simplify(out);
        end
        function out=exp(qN,X,n,inv)
            arguments
                qN
                X
                n
                inv (1,1) logical=false
            end
            q=qN.q;
            if inv
                q=q^-1;
            end
            if nargin==2
                n=5;
            end
            out=0;
            for k=0:n
                out=out+(q^(k*(k-1)/2)/qN.fac(k))*X^k;
            end
            % out=1;
            % for k=1:n
            %     out=X*(q^(n-k)/qN.n(n+1-k))*out+1;
            % end
        end
        function ret=binomial(qN,n,m)
            % out=simplify(qN.fac(n)/(qN.fac(n-m)*qN.fac(m)));
            if ~isa(n,'sym')&&n<m
                ret=0;
                return
            else
                ret=1;
            end
            for ii=0:m-1
                ret=qN.n(n-ii)/qN.n(m-ii)*ret;
            end
        end
        function ret=factor(qN,arg)
            F=factor(sym(arg));
            N=numel(F);
            T=table(nan(N,1),F.',VariableNames={'type','value'});
            ret=T;
            q=qN.q;
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


    end

end