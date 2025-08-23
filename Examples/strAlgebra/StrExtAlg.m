classdef(InferiorClasses=?sym) StrExtAlg<StrAlg
    properties(Constant,Hidden)
        B=TypeParam(@(N)Bases(N,"X"+(1:N),"Ext"+N))
    end
    properties
        Nvar
    end

    %% generation
    methods
        function obj=StrExtAlg(Nvar)
            obj.Nvar=Nvar;
            obj.base=obj.B.get(Nvar);
        end
        % function ret=getActMono(obj,idxTesor,idxFactor)
        %     pw=obj.pw{idxTesor}(idxFactor);
        %     N=obj.Nvar;
        %     if any(1:N==pw)
        %         obj=obj.make("x",pw);
        %     elseif any(1:N==pw-N)
        %         obj=obj.make("x",pw-N,-1);
        %     elseif any(1:N==pw-2*N)
        %         obj=obj.make("qth",pw-2*N,1,obj.q);
        %     elseif any(1:N==pw-3*N)
        %         obj=obj.make("qth",pw-3*N,-1,obj.q);
        %     end
        %     ret=obj.getActMono@StrEndV(1,1);
        % end
        function obj=unit(obj)
            obj=unit@StrAlg(obj);
        end
        function obj=make(obj,arg,varargin)
            obj=obj.make@StrAlg(arg,varargin{1},obj.base);
        end
        function ret=mpower(i1,i2)
            if isa(i2,'double')&&i2<0
                assert(i1.term==1,"minus power not allowed")
                i2=-i2;
                i1=1/i1;
            end
            ret=mpower@StrAlg(i1,i2);
        end
        function ret=mrdivide(i1,i2)
            assert(isa(i2,"StrWeylXQ"))
            assert(i2.term==1,'除算が定義されません')
            Nvar=i2.Nvar;
            arr=(1:Nvar)'+[1 0 3 2]*Nvar;
            arr=arr(:)';
            i2inv=i2.set_cp(i2.cf,cellfun(@(p){arr(flip(p))},i2.pw),i2.bs);
            ret=i1*i2inv;
        end
        function ret=Delta(obj)
            I=obj.unit;
            ret=obj.algfun(@fun,I|I);
            function ret=fun(p,b)
                ret=(I.make(1,{p})|I)+(I|I.make(1,{p}));
            end
        end
    end
    methods(Static)
        function [O,X]=getGenerator(Nvar)
            arguments
                Nvar
            end
            O=StrExtAlg(Nvar);
            O=O.make(0,{[]});
            X=dictionary();
            for i=1:Nvar
                X(i)=O.make(1,{i});
            end
        end
    end
    methods
        function ret=convertBaseString(obj,pw,bs)
            ret=convertBaseString@StrAlg(obj,pw,bs);
        end
        %% relation
        function ret=replace(obj,Ntimes)
            idx=all(cellfun(@allunique,obj.pw),2);
            obj.cf=obj.cf(idx);
            obj.pw=obj.pw(idx,:);
            ret=replace@StrAlg(obj,Ntimes);
        end
        function ret=mtimes(i1,i2)
            [i1,i2]=alignNum(i1,i2);
            assert(isequal(i1.base,i2.base),'異なる空間での積エラー')
            R=i1.rank;
            ret=lfun_(i1|i2,@fun);
            ret.base=i1.base;
            ret.ZERO=i1.ZERO;
            ret=ret.calc();
            function [c,p,b]=fun(p,b)
                c=1;
                lp=cellfun(@length,p);
                for j=1:R
                    for k=1:j-1
                        c=c*(-1)^(lp(j)*lp(k+R));
                    end
                end
                p=cellfun(@(p1,p2)[p1 p2],p(1:R),p(R+1:end),UniformOutput=false);
            end
        end
        function [rel,mlist,comm,inv]=get2vRelation(arg)
            persistent RS
            if isempty(RS)
                RS=TypeParam(@createRel);
            end
            S=RS.get(arg.Nvar);
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
            function S_=createRel(Nvar)
                O=StrExtAlg.getGenerator(Nvar);
                S_=struct;
                % xy+yx=0
                S_.rel=O.empty();

                for i=2:Nvar
                    for j=1:i-1
                        S_.rel(end+1)=O.make([1 1],{[i j],[j i]});
                    end
                end
                % x^2=0はparentで簡約化しない
                S_.comm=[nan;nan];
                S_.inv=[nan;nan];
                S_=O.get2vRelation_(S_);
                function ret=expr1(k1,k2)
                    % 反交換関係
                    ret=O.make([1 1],{[k1,k2],[k2,k1]});
                end
            end
        end

    end
end
function ret=i2p(arg)
    ret=2.^(arg-1)+1;
end