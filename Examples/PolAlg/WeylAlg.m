classdef WeylAlg<PolAlg&IAlg
    %WEYLALG このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        Nvar=0
    end

    methods
        function out=mtimes(i1,i2)
            arguments
                i1 WeylAlg
                i2 WeylAlg
            end
            NvarMax=max(i1.Nvar,i2.Nvar);
            if i1.Nvar<NvarMax
                i1=i1.setBase(i2.base);
            elseif i2.Nvar<NvarMax
                i2=i2.setBase(i1.base);
            end
            base=i2.base;
            out2=i1|i2;
            for i0=1:NvarMax
                out2=out2.lfun(@multiply);
            end
            out=out2.prodeval;
            out.base=base;
            function [c,p]=multiply(p)
                i=i0;
                n=NvarMax;
                Kmax=min(p([n+i 2*n+i]));
                K=(0:Kmax)';
                c=arrayfun(@(k)nchoosek(p(n+i),k)*nchoosek(p(2*n+i),k)*factorial(k),K);
                % p=(p(1:2)+p(3:4))-K;
                p=repmat(p,Kmax+1,1);
                p(:,[i n+i ])=p(:,[i n+i])+p(:,[2*n+i 3*n+i])-K;
                p(:,[2*n+i 3*n+i])=0;
            end
        end
        function obj=setBase(obj,base)
            Nvar=base.dim/2;
            if obj.Nvar==Nvar
                return
            elseif obj.Nvar==0
                obj.base=base;
                % obj.pw=[obj.pw obj.ptype.zeros(obj.term,2*(Nvar-obj.Nvar))];
            else
                error down_casting
            end
            obj.Nvar=Nvar;
            % obj.base=Bases(2*Nvar,["X"+(1:Nvar),"D"+(1:Nvar)]);
        end


    end
    methods(Static)
        function [O,X,D]=getGenerator(Nvar,name)
            if ~exist("name","var")
                name=["X"+(1:Nvar),"D"+(1:Nvar)];
            end
            O=WeylAlg(0).setBase(Bases(2*Nvar,name));

            X=repmat(O,1,Nvar);
            D=repmat(O,1,Nvar);
            for i=1:Nvar
                X(i).pw(1,i)=1;
                D(i).pw(1,i+Nvar)=1;
            end
            [X.cf,D.cf]=deal(1);
        end


    end
end

