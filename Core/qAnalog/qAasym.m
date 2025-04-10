classdef qAasym
    properties
        q
    end
    methods

        function out=qAasym(q)
            out.q=q;
        end
        function out=eps(in,q)
            if ~exist("q","var")
                q=in.q;
            end
            out=1-q;
        end
        function out=num(in,n,q)
            if ~exist("q","var")
                q=in.q;
            end
            out=(1-q^n)/(1-q);
        end
        function out=pch(in,z,n,q)
            if ~exist("q","var")
                q=in.q;
            end
            if ~exist("n","var")
                n=10;
            end
            outs=(in.q-in.q+1)*ones(1,numel(z));
            for ii=1:numel(outs)
                for kk=0:n-1
                    outs(ii)=outs(ii)*(1-z(ii)*q^kk);
                end
            end
            out=prod(outs);
        end
        function out=fac(in,n,q)
            if ~exist("q","var")
                q=in.q;
            end
            out=in.pch(q,n,q)/in.eps(q)^n;
        end
        function out=exp(in,z,n,q)
            if ~exist("q","var")
                q=in.q;
            end
            if ~exist("n","var")
                n=10;
            end
            out=in.q-in.q;
            for kk=0:n
                out=out+1/in.fac(kk,q)*z^kk;
            end
        end
        function out=phi(in,a,b,z,n,q)
            if ~exist("q","var")
                q=in.q;
            end
            if ~exist("n","var")
                n=5;
            end
            modify=1+numel(b)-numel(a);
            out=in.q-in.q;
            for kk=0:n
                out=out+in.pch(a,kk,q)/in.pch([b(:);q],kk,q)...
                    *((-1)^kk*q^(kk*(kk-1)/2))^modify*z^kk;
            end

        end
    end
end