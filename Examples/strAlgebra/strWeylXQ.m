classdef(InferiorClasses=?sym) StrWeylXQ<StrEndV
    properties(Constant,Hidden)
        B=TypeParam(@(N)Bases(4*N,reshape(["X","Xi","qt","qti"]+(1:N)',1,4*N),"weyl_"+N))
    end
    properties
        Nvar
        q=sym('q')
    end

    %% generation
    methods
        function obj=StrWeylXQ(Nvar)
            obj.Nvar=Nvar;
            obj.algbase=obj.B.get(Nvar);
        end
        function ret=getActMono(obj,idxTesor,idxFactor)
            pw=obj.pw{idxTesor}(idxFactor);
            N=obj.Nvar;
            if any(1:N==pw)
                obj=obj.make("x",pw);
            elseif any(1:N==pw-N)
                obj=obj.make("x",pw-N,-1);
            elseif any(1:N==pw-2*N)
                obj=obj.make("qth",pw-2*N,1,obj.q);
            elseif any(1:N==pw-3*N)
                obj=obj.make("qth",pw-3*N,-1,obj.q);
            end
            ret=obj.getActMono@StrEndV(1,1);
        end
        function obj=unit(obj)
            obj=unit@StrAlg(obj);
        end
        function obj=make(obj,arg,varargin)
            if isString(arg)
                obj=obj.make@StrEndV(arg,varargin{:});
            else
                obj=obj.make@StrAlg(arg,varargin{1},obj.algbase);
            end
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
    end
    methods(Static)
        function [O,X,qth]=getGenerator(Nvar)
            arguments
                Nvar
            end
            O=StrWeylXQ(Nvar);
            O=O.make(0,{[]});
            C=num2cell(1:Nvar);
            C(2,:)=cellfun(@(i){O.make(1,{i})},C(1,:));
            C(3,:)=cellfun(@(i){O.make(1,{2*Nvar+i})},C(1,:));
            X=dictionary(C{[1 2],:});
            qth=dictionary(C{[1 3],:});

        end
    end
    methods
        function ret=convertBaseString(obj,pw,bs)
            ret=convertBaseString@StrAlg(obj,pw,bs);
        end
        %% relation
        function ret=replace(obj,Ntimes)
            ret=replace@StrAlg(obj,Ntimes);
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
                O=StrWeylXQ.getGenerator(Nvar);
                S_=Struct;
                % q^θ X=qXq^θ
                S_.rel=[expr1(0,2),expr1(2,1),expr1(3,0),expr1(1,3)];
                        % expr2(0,1),expr2(1,0),expr2(2,3),expr2(3,2)];
                
                T=combinations(1:4*Nvar,1:4*Nvar);
                T2=[expr3(0,2),expr3(1,2),expr3(0,3),expr3(1,3)];
                T3=setdiff(T{:,:},T2','rows');
                S_.comm=T3(T3(:,1)>T3(:,2),:)';
                S_.inv=[0;Nvar]+[1:Nvar,2*Nvar+(1:Nvar)];
                S_=O.get2vRelation_(S_);
                function ret=expr1(k1,k2)
                    % q可換な関係式を生成 q^θ x=q x q^θ
                    ret=arrayfun(@(i)O.make([1 -arg.q], ...
                        {[i+Nvar*k2 i+Nvar*k1] [i+Nvar*k1 i+Nvar*k2]}) ...
                        ,1:Nvar);
                end
                function ret=expr2(k1,k2)
                    % 逆元の関係式 x x^-1=1
                    ret=arrayfun(@(i)O.make([1 -1],{i+Nvar*[k1 k2] []}),1:Nvar);
                end
                function ret=expr3(k1,k2)
                    % 逆元の関係式 x x^-1=1
                    ret=[k2;k1]*Nvar+(1:Nvar);
                end
            end
        end

    end
end