classdef(InferiorClasses=?sym) VectExtAlg<VectAlg
    % uqsl2borelsmall borel subalgebra of Uq(sl2)
    % style=R: Δ(E)=1⊗E+E⊗K, style=L: Δ(E)=E⊗1+K⊗E
    properties(Constant)
        bs0= TypeParam(@makeBase)
    end
    properties
        n
    end
    methods(Static)
        function [Z,X]=getGenerator(n)
            Z=VectExtAlg();
            Z.n=n;
            Z=Z.setBase(Z.bs0.get(n));
            Z=Z.make([],[]);

            
            Z.setConst();
            X=dictionary;
            for i=1:Z.n
                X(i)=Z.make(1,2^(i-1)+1);
            end
        end
    end
    methods
        function obj=simplify(obj)
            cf=obj.cf;
            P=obj.p;
            for i=1:numel(cf)
                [N,D]=rat(cf(i));
                cf(i)=mod(N*D^(P-2),P);
            end
            obj.cf=cf;
        end
        function obj=casttype(obj,arg)
            if isa(arg,'double')&&isscalar(arg)
                obj=obj.set_c([arg;zeros(obj.p-1,1)]);
            else
                error not_impl
            end
        end
        function setConst(obj)
            N=3;
            Z=StrExtAlg.getGenerator(N);
            b=fliplr(dec2bin(0:2^N-1)=='1');
            powerset = arrayfun(@(k){find(b(k,:))}, 1:2^N);
            powersetC=cellfun(@char,powerset,Un=0);
            M = zeros(2^N*[1,1,1]);
            for i=1:2^N
                for j=1:2^N
                    P=calc(Z.make(1,powerset(i))*Z.make(1,powerset(j)));
                    cf=P.cf;
                    if cf==0,continue; end
                    pw={char(P.pw{1})};
                    [~,idx]=ismember(pw,powersetC);
                    M(i,j,idx)=cf;
                end
            end
            C=zeros(2^N*[1,1,1]);
            for i=1:2^N
                D=calc(Delta(Z.make(1,powerset(i))));
                cf=D.cf;
                pw=cellfun(@char,D.pw,Un=0);
                [~,idx]=ismember(pw,powersetC);
                for j=1:size(idx,1)
                    C(i,idx(j,1),idx(j,2))=cf(j);
                end
            end
            epsilon=zeros(2^N,1);1;
            eta = epsilon;
            lp=cellfun(@length,powerset);
            S=diag((-1).^lp);
            obj.setSC(obj.identifier,M,eta,C,epsilon,S);
        end
    end
end
function ret=makeBase(n)
    name=['ExtAlg' num2str(n)];
    b=fliplr(dec2bin(0:2^n-1)=='1');
    xs="x"+(1:n);
    var=strings(1,2^n);
    for i=2:2^n
        var(i)=join(xs(b(i,:)),"");
    end
    var(1)="1";
    ret=Bases(2^n,var,name);
    
end