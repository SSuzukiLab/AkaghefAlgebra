classdef QWeylAlg<symp&IAlg
    %WEYLALG このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        Nvar=0
        q
    end

    methods
        function obj=QWeylAlg(varargin)
            obj@symp(varargin{:})
            obj.ctype="S";
        end
        function arg=prodevalQW1(arg)
            % 新しい公式で計算
            x=arg;
            n=1;
            q=arg.q;
            for i=1:n
                arg=arg.lfun(@multiply);
            end

            function [c,p]=multiply(p)
                pw_=p([n+i 2*n+i]);
                pws=sum(pw_);
                k=min(pw_);
                qN=qNumS(q);
                nt=(k:-1:0)';
                pw2=pw_+0*nt;
                pw_=pw_-nt;
                c=prod(arrayfun(@qN.fac,pw2) ...
                    ./arrayfun(@qN.fac,pw_),2)...
                    ./arrayfun(@qN.fac,nt)...
                    .*q.^(-nt.*(nt-1)/2)...
                    .*q.^(2*prod(pw2,2))...
                    .*q.^(-pws*nt+nt.*(nt-1)/2*2);
                % ret={simplify(expr.cf./cf).'};
                p=pw_+[p([i 3*n+i])];
            end
        end
        function arg=prodevalQW2(arg)
            % ワイル代数の関係式から逐次計算
            arguments
                arg QWeylAlg
            end
            x=arg;
            n=1;
            for i0=1:n
                arg=arg.lfun(@multiply);
            end
            function [c,p]=multiply(p)
                i=i0;
                p2=[p(n+i),p(2*n+i)];
                if min(p2)==0
                    c=1+arg.ctype.zero;
                    p=p(1:2*n)+p(2*n+1:4*n);
                else
                    [dif1,dif2]=deal(zeros(1,4*n));
                    dif1([n+i 3*n+i])=[-1 +1];
                    dif2([n+i,2*n+i])=[-1 -1];
                    px=p(2*n+i);
                    x=x.set_cp([x.q^(2*px);x.q^(px-1)*qNumS(x.q).n(px)],[p+dif1;p+dif2]);
                    x=x.prodevalQW2;
                    c=x.cf;
                    p=x.pw;
                end
            end
        end

        function out=mtimesQW1(i1,i2)
            % 正規順序積の明示公式により計算
            arguments
                i1 QWeylAlg
                i2 QWeylAlg
            end
            NvarMax=max(i1.Nvar,i2.Nvar);
            if i1.Nvar<NvarMax
                i1=i1.setBase(i2.base);
            elseif i2.Nvar<NvarMax
                i2=i2.setBase(i1.base);
            end
            out2=i1|i2;
            out=out2.prodevalQW1;
            out.base=i2.base;
        end
        function out=mtimesQW2(i1,i2)
            arguments
                i1 QWeylAlg
                i2 QWeylAlg
            end
            NvarMax=max(i1.Nvar,i2.Nvar);
            if i1.Nvar<NvarMax
                i1=i1.setBase(i2.base);
            elseif i2.Nvar<NvarMax
                i2=i2.setBase(i1.base);
            end
            base=i2.base;
            out2=i1|i2;
            out=out2.prodevalQW2;
            out.base=i2.base;
        end
        function out=mtimes(i1,i2)
            out=mtimesQW1(i1,i2);
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
        function [O,X,D]=getGenerator(Nvar,name,q)
            arguments
                Nvar
                name =["X"+(1:Nvar),"D"+(1:Nvar)]
                q =sym('q')
            end
            assert(Nvar==1,"un implemented")
            O=QWeylAlg(0).setBase(Bases(2*Nvar,name));
            O.q=q;
            O.ctype="S";
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

