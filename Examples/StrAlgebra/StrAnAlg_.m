classdef(InferiorClasses=?sym) StrAnAlg_<StrAlg&HopfAlg
    properties(Constant,Hidden)
        B=TypeParam(@(N)Bases(4*N,reshape(["E","F","K","Ki"]+(1:N)',1,4*N),"Uqsl_"+(N+1)))
    end
    properties
        Nvar
        q=sym('q')
    end

    %% generation
    methods
        % function obj=StrAnalg(Nvar)
        %     obj.Nvar=Nvar;
        %     obj.algbase=obj.B.get(Nvar);
        % end
        % function obj=unit(obj)
        %     obj=unit@StrAlg(obj);
        % end
        function obj=make(obj,cf,pw)
            obj=obj.make@StrAlg(sym(cf),pw,obj.algbase);
        end
    end
    methods(Static)
        function [O,E,F,K,Ki]=getGenerator(Nvar)
            arguments
                Nvar {mustBeInteger}
            end
            O=StrAnAlg_();
            O.ctype="S";
            O.Nvar=Nvar;
            O.algbase=O.B.get(Nvar);
            O=O.make(0,{[]});
            C=num2cell(1:Nvar);
            C(2:5,:)=reshape(arrayfun(@(i){O.make(1,{i})},1:4*Nvar),Nvar,4).';
            E=dictionary(C{[1 2],:});
            F=dictionary(C{[1 3],:});
            K=dictionary(C{[1 4],:});
            Ki=dictionary(C{[1 5],:});


        end
    end
    methods
        function ret=convertBaseString(obj,pw,bs)
            ret=convertBaseString@StrAlg(obj,pw,bs);
        end
        %% relation
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
                mustBeInteger(Nvar)
                O=StrAnAlg_.getGenerator(Nvar);
                % S=struct;
                q=sym('q');
                S_=struct;
                % q^θ X=qXq^θ
                % [K,E], [Ki,E] KKi [E,F]=[H],[Ei,Fj]=0 [K,Ki]
                CM=diag(2*ones(1,Nvar))-diag(ones(1,Nvar-1),1)-diag(ones(1,Nvar-1),1);
                for i=1:Nvar
                    for j=1:Nvar
                        if CM(i,j)==-1

                        end
                    end
                end
                T=combinations(1:Nvar,1:Nvar);
                T2=T;
                T2(T2.Var1==T2.Var2,:)=[];
                Ti=abs(T2.Var1-T2.Var2)==1;
                T3=T2(~Ti,:);
                T4=T2(Ti,:);
                T5=table((1:Nvar)',(1:Nvar)');
                S_.rel=[expr1(T5,0,2,2),expr1(T5,1,2,-2),expr1(T5,0,3,-2),expr1(T5,1,3,2), ...
                    expr1(T4,0,2,-1),expr1(T4,1,2,1),expr1(T4,0,3,1),expr1(T4,1,3,-1), ...
                    expr2(2,3),expr3(),expr4(T4),expr4(T4+Nvar)];
                S_.comm=[[T2.Var1,T2.Var2+Nvar];
                    [T.Var1+Nvar*2,T.Var2+Nvar*2];
                    [T.Var1+Nvar*3,T.Var2+Nvar*3];
                    [T.Var1+Nvar*2,T.Var2+Nvar*3];
                    [T.Var1+Nvar*3,T.Var2+Nvar*2];
                    [T3.Var1+Nvar*0,T3.Var2+Nvar*2];
                    [T3.Var1+Nvar*0,T3.Var2+Nvar*3];
                    [T3.Var1+Nvar*1,T3.Var2+Nvar*2];
                    [T3.Var1+Nvar*1,T3.Var2+Nvar*3];
                    ]';
                S_.inv=[nan;nan];
                S_=O.get2vRelation_(S_);
                function ret=expr1(T,k1,k2,n)
                    % q可換な関係式を生成 KE=q^nEK,KF=q^-nFK
                    k1=Nvar*k1;
                    k2=Nvar*k2;
                    ret=arrayfun(@(i,j)O.make([1 -q^n], ...
                        {[k2+j k1+i] [k1+i k2+j]}) ...
                        ,T{:,1}',T{:,2}');
                end
                function ret=expr2(k1,k2)
                    % 逆元の関係式 x x^-1=1
                    ret=arrayfun(@(i)O.make([1 -1],{i+Nvar*[k1 k2] []}),1:Nvar);
                end
                function ret=expr3()
                    % 関係式を生成 [E,F]=[H]
                    ret=arrayfun(@(i)O.make([1 -1 -[q -q]/(q^2-1)], ...
                        {[i+Nvar*0 i+Nvar*1] [i+Nvar*1 i+Nvar*0] i+Nvar*2 i+Nvar*3}) ...
                        ,1:Nvar);
                end
                function ret=expr4(T)
                    % Serre関係式を生成 ad(E)^2=0,ad(E)^2=0
                    ret=arrayfun(@(i,j)O.make([1 -(q+q^-1) 1], ...
                        {[i i j] [i j i] [j i i]}) ...
                        ,T{:,1}',T{:,2}');
                end
            end
        end
        %% rep

        function ret = repMono(obj)
            persistent RS
            if isempty(RS)
                RS=TypeParam(@createRep);
            end
            arr=RS.get(obj.Nvar);
            I=arr{1}(end);
            converted=obj.algfun(@fun,I);
            ret=I.set_cp(converted.cf,converted.pw,converted.bs);
            ret.dimV=obj.Nvar+1;
            function ret=fun(p,b)
                % assert(b==StrStrUqsl2.B)
                ret=arr{1}(p);
            end
            function ret=createRep(N)
                q=obj.q;
                [O,x,qth]=StrWeylXQ.getGenerator(N+1);
                I_=O.unit;
                arr_=[arrayfun(@(k)(x(k)/x(k+1))*(qth(k+1)-1/qth(k+1))*(q-q^-1)^-1,1:N), ...
                    arrayfun(@(k)(x(k+1)/x(k))*(qth(k)-1/qth(k))*(q-q^-1)^-1,1:N), ...
                    arrayfun(@(k)qth(k)/qth(k+1),1:N), ...
                    arrayfun(@(k)qth(k+1)/qth(k),1:N), ...
                    I_];
                ret={arr_};
            end
        end
        function v=getModGenerator(obj)
            N=obj(1).Nvar+1;
            v0=PolAlg(1,zeros(1,N));
            v0.base.ctype="S";
            v0.ctype="S";
            v=cellfun(@(p)v0.set_cp(1,p),mat2cell(eye(N),ones(1,N),N));
        end
        %% hopf
        function ret = Delta(obj)
            persistent List
            if isempty(List)
               List=TypeParam(@createDelta);
            end
            arr=List.get(obj.Nvar);
            I=arr{1}(end);
            ret=obj.algfun(@fun,I|I);
            function ret=fun(p,b)
                % assert(b==StrStrUqsl2.B)
                ret=arr{1}(p);
            end
            function ret=createDelta(N)
                 [O,E,F,K,Ki]=obj.getGenerator(N+1);
                I_=O.unit;
                arr=[arrayfun(@(n)(E(n)|I_)+(K(n)|E(n)),1:N),...
                    arrayfun(@(n)(F(n)|Ki(n))+(I_|F(n)),1:N),...
                    arrayfun(@(n)K(n)|K(n),1:N),...
                    arrayfun(@(n)Ki(n)|Ki(n),1:N),...
                    I_];
                ret={arr};
            end
        end

        % Counit
        function ret = counit(obj)
            ret=obj.lfun(@fun);
            ret=sum(ret.cf);
            function [c,p,b]=fun(p,b)
                c=all(p{1}==0);
                p={[]};
                b={[]};
            end
        end

        % Antipode
        function ret = antipode(obj)
            % NG
            ret=obj.lfun(@fun);
            function [c,p,b]=fun(p,b)
                c=(-1)^length(p{1});
            end
        end
    end
end