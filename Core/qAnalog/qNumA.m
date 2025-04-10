classdef qNumA
    properties
        q
    end
    methods

        function out=qNumA(q)
            out.q=q;
        end
        function out=eps(in)
            out=1-in.q^-2;
        end
        function out=num(in,n)
            out=(in.q^n-1)/(in.q-1);
        end
        function out=pch(in,z,n)
            outs=ones(1,numel(z),class(in.q));
            for ii=1:numel(outs)
                for kk=0:n-1
                    outs(ii)=outs(ii)*in.num(z(ii)+kk);
                end
            end
            out=prod(outs);
        end
        function out=fac(in,n)
            out=in.pch(1,n);
        end
        function out=exp(in,z,n)
            q=in.q;
            out=q-q;
            for kk=0:n
                out=out+(q^(kk*(kk-1)/2)/in.fac(kk,q))*z^kk;
            end
        end
        function ret=nchoosek(in,a,b)
            q=in.q;
            ret=q-q+1;
            for i=0:b-1
                ret=ret*((1-q^(a-i))/(1-q^(b-i)));
            end
        end
        function out=phi(in,a,b,z,n)
            val=numel(a)-numel(b)-1;
            Z=z/q^(sum(a)-sum(b)-1);
            out=in.q-in.q;
            for kk=0:n
                out=out+in.pch(a,kk,q)/in.pch([b(:);1],kk,q)...
                    *((q^-1-q)^kk*q^(kk*(kk-1)/2))^val*Z^kk;
            end
        end
    end
end