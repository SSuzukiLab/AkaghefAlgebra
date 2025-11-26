classdef qNumA
    methods(Static)
        
        % q-epsilon: out = 1 - q^-2
        function out=eps(q)
            out=1-q^-2;
        end
        
        % q-number [n]_q: out = (q^n - 1) / (q - 1)
        function out=n(q,n,is_series)
            if nargin<3||~is_series
                out=(q^n-1)/(q-1);
            else
                out=0;
                for kk=0:n-1
                    out=out+q^kk;
                end
            end     
        end
        
        % q-Pochhammer symbol (z; q)_n
        % オリジナル: in.num(z(ii)+kk) の使用を確認
        function out=pch(q,z,n)
            out=ones(size(z),class(q));
            assert(isscalar(n))
            if n<0
                out=1./qNumA.pch(q^-1,z*q,-n);
                return
            end
            for ii=1:numel(out)
                out(ii)=prod(1-z(ii)*q.^(0:n-1));
            end
        end
        
        % q-factorial [n]!_q: out = (1; q)_n
        % pch(q, 1, n) を使用
        function out=fac(q,n)
            out=qNumA.pch(q,1,n);
        end
        
        % q-exponential E_q(z)
        % pch と fac が q を第一引数として取るよう変更したため、呼び出しも修正
        function out=exp(q,z,n)
            out=q-q; % 0を生成
            for kk=0:n
                % fac(q, kk) を使用
                out=out+(q^(kk*(kk-1)/2)/qNumA.fac(q,kk))*z^kk; 
            end
        end
        
        % q-binomial coefficient [a choose b]_q
        function ret=binom(q,a,b)
            ret=q-q+1; % 1を生成
            for i=0:b-1
                ret=ret*((1-q^(a-i))/(1-q^(b-i)));
            end
        end
        
        % q-hypergeometric function: _r phi _s (a; b; q, z)_n
        % pch が q を第一引数として取るよう変更したため、呼び出しも修正
        function out=phi(q,a,b,z,n)
            val=numel(a)-numel(b)-1;
            Z=z/q^(sum(a)-sum(b)-1);
            out=q-q; % 0を生成
            for kk=0:n
                % pch(q, a, kk) と pch(q, [b(:);1], kk) を使用
                out=out+qNumA.pch(q,a,kk)/qNumA.pch(q,[b(:);1],kk)...
                    *((q^-1-q)^kk*q^(kk*(kk-1)/2))^val*Z^kk;
            end
        end
    end
end