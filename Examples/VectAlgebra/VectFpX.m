classdef(InferiorClasses=?sym) VectFpX<VectAlg
    % uqsl2borelsmall borel subalgebra of Uq(sl2)
    % style=R: Δ(E)=1⊗E+E⊗K, style=L: Δ(E)=E⊗1+K⊗E
    properties(Constant)
        bs0= TypeParam(@makeBase)
    end
    properties
        p
    end
    methods(Static)
        function [Z,X]=getGenerator(p)
            Z=VectFpX();
            Z=Z.setBase(Z.bs0.get(p));
            Z.p=p;
            Z.setConst();
            X=Z.make(1,2);
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
            P=obj.p;
            M = zeros(P*[1 1 1]);
            for i=0:P-1
                for j=0:P-1
                    if i+j>=P
                        continue
                    end
                    M(i+1,j+1,i+j+1) =1;
                end
            end
            C = zeros(P*[1 1 1]);
            for i=0:P-1
                for k=0:i
                    C(i+1,k+1,i-k+1) = nchoosek(i,k);
                end
            end
            [epsilon,eta] =deal(zeros(P,1));
            epsilon(1)=1;
            eta(1)=1;
            eta = [1; zeros(P-1,1)];
            S=diag((-1).^(0:P-1));
            obj.setSC(obj.identifier,M,eta,C,epsilon,S);
        end
    end
end
function ret=makeBase(p)
    name=['FpX' num2str(p)];
    var="x^"+(0:p-1);
    var(1)="1";
    ret=Bases(p,var,name);

end