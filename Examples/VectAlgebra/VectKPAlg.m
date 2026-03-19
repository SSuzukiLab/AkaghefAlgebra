classdef(InferiorClasses=?sym) VectKPAlg<VectAlg
    % Kac-Paljutkin Algebra
    % https://arxiv.org/pdf/2505.00645
    properties(Constant)
        bs0=makeBase()
    end
    methods(Static)
        function [Z,x,y,z]=getGenerator()
            Z=VectKPAlg();
            Z=Z.setBase(Z.bs0);
            Z.setConst;
            x=arrayfun(@(p)Z.make(1,p),[2,3,5]);
            if nargout>2
                [x,y,z]=deal(x(1),x(2),x(3));
            end
        end
    end
    methods
        function obj=casttype(obj,arg)
            if isa(arg,'double')&&isscalar(arg)
                obj=obj.set_c([arg;zeros(7,1)]);
            else
                error not_impl
            end
        end
        function [M,C,eta,epsilon,S] = setConst(obj)
            %KP_H8_structure_constants  Structure constants of Kac–Paljutkin Hopf algebra H8.
            O=StrKPAlg.getGenerator();
            % P:powers of algebra generators
            % X:basis of the algebra
            P=combinations([0,1],[0,1],[0,1]);
            P=fliplr(P{:,:});
            for i=1:size(P,1)
                tmp=arrayfun(@(i,n){repmat(i,n)},1:size(P,2),P(i,:));
                pw(i)={horzcat(tmp{:})};
                X(i)=O.make(1,pw(i));
            end
            eidx=find(cellfun(@(p)isempty(p),pw));
            if isscalar(eidx)
                pw([eidx,end+1,end+2])={[],zeros(1,0),zeros(0,1)};
            end
            D=length(X);
            M=zeros(D,D,D);
            eta=zeros(D,1);
            C=zeros(D,D,D);
            epsilon=zeros(D,1);
            S=zeros(D,D);
            S2V=dictionary(pw,[1:D,eidx,eidx]);
            tmp=O.unit;
            eta(S2V(tmp.pw))=tmp.cf;
            for i=1:D
                for j=1:D
                    tmp=X(i)*X(j);
                    M(i,j,S2V(tmp.pw))=tmp.cf;
                end
                tmp=Delta(X(i));
                for j=1:tmp.term
                    C(i, S2V(tmp.pw(j,1)),S2V(tmp.pw(j,2)) )=tmp.cf(j);
                end
                tmp=counit(X(i));
                epsilon(i)=tmp; % tmp.cf?
                tmp=antipode(X(i));
                S(i,S2V(tmp.pw))=tmp.cf;
            end
            obj.setSC(obj.identifier,M,eta,C,epsilon,S);
        end

    end
end
function b=makeBase()
    bstr=strrep(string(char((fliplr(dec2bin(0:7))=='1').*('x'+(0:2)))),char(0),'').';
    bstr(1)="1";
    b=Bases(8,bstr,"KP8");
end