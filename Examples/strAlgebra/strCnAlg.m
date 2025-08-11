classdef(InferiorClasses=?sym) strCnAlg<strQUEAlg
    properties(Constant,Hidden)
        % B
        B            =TypeParam(@(N)Bases(5*N,reshape(["E","F","K","Ki" "H"]+(1:N)',1,5*N),"Uqsl_"+(N+1)))
        CM           =TypeParam(@createDelta)
        RVL          =TypeParam(@createRVL)
        DeltaStorage =TypeParam(@createDelta)
        RepStorage   =TypeParam(@createRep)
        RelStorage   =TypeParam(@createRel)
        Vtype='symCnmodalg'
    end
    properties(Dependent)
        dimV
        
    end

    %% generation
    methods
        % function obj=strAnalg(Nvar)
        %     obj.Nvar=Nvar;
        %     obj.algbase=obj.B.get(Nvar);
        % end
        % function obj=unit(obj)
        %     obj=unit@strAlg(obj);
        % end
        function obj=make(obj,cf,pw)
            obj=obj.make@strAlg(sym(cf),pw,obj.algbase);
        end
        function ret=get.dimV(obj)
            ret=obj.Nvar+1;
        end
    end
    methods(Static)
        function [O,E,F,K,Ki]=getGenerator(Nvar)
            arguments
                Nvar {mustBeInteger}
            end
            O=strCnAlg();
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
            ret=convertBaseString@strAlg(obj,pw,bs);
        end
        %% relation

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
            ret.dimV=obj.Nvar*2;


            function ret=fun(p,b)
                % assert(b==Uqsl2.B)
                ret=arr{1}(p);
            end
            function ret=createRep(N)
                q=obj.q;
                d=dictionary([-N:-1 1:N],1:2*N);
                [O,x,qth]=strWeylXQ.getGenerator(N*2);
                I_=O.unit;
                x=dictionary(d.keys,x(d.values));
                qth=dictionary(d.keys,qth(d.values));
                Y=repmat(O,[N 5]);
                Y(1,1)=(x(-1)/x(1))*(qth(1)^2-1/qth(1)^2)*(q^2-q^-2)^-1;
                Y(1,2)=(x(1)/x(-1))*(qth(-1)^2-1/qth(-1)^2)*(q^2-q^-2)^-1;
                for j=2:N
                    q_expr=I_;
                    for k=1:j-1
                        q_expr=q_expr*qth(k)*qth(-k);
                    end

                    Y(1,1)=Y(1,1)+((q^2-1)/(q^2+1))*(x(-j)*x(j)/x(1)^2)*...
                    (q^-1*qth(1)-q/qth(1))*(qth(1)-1/qth(1))*(q-q^-1)^-2*q_expr;
                    Y(1,2)=Y(1,2)+((q^2-1)/(q^2+1))*(x(-j)*x(j)/x(-1)^2)*...
                    (q^-1*qth(-1)-q/qth(-1))*(qth(-1)-1/qth(-1))*(q-q^-1)^-2*q_expr;
                end
                Y(1,3)=qth(-1)^2/qth(1)^2;
                Y(1,4)=qth(1)^2/qth(-1)^2;
                for i=2:N
                    Y(i,1)=(x(-i)/x(1-i))*(qth(1-i)-1/qth(1-i))*(q-q^-1)^-1*qth(i-1)/qth(i) ...
                        -(x(i-1)/x(i))*(qth(i)-1/qth(i))*(q-q^-1)^-1;
                    Y(i,2)=(x(1-i)/x(-i))*(qth(-i)-1/qth(-i))*(q-q^-1)^-1 ...
                        -(x(i)/x(i-1))*(qth(i-1)-1/qth(i-1))*(q-q^-1)^-1*qth(1-i)/qth(-i);
                    Y(i,3)=qth(-i)/qth(1-i)*qth(i-1)/qth(i);
                    Y(i,4)=1/qth(-i)*qth(1-i)/qth(i-1)*qth(i);
                end

                % Y=[arrayfun(@(k)(x(k)/x(k+1))*(qth(k+1)-1/qth(k+1))*(q-q^-1)^-1,1:N), ...
                %     arrayfun(@(k)(x(k+1)/x(k))*(qth(k)-1/qth(k))*(q-q^-1)^-1,1:N), ...
                %     arrayfun(@(k)qth(k)/qth(k+1),1:N), ...
                %     arrayfun(@(k)qth(k+1)/qth(k),1:N), ...
                %     I_];
                ret={[Y(:).' I_]};
            end


            function v=getModGenerator(obj)
                N=obj(1).Nvar;
                v0=symp(1,zeros(1,N));
                v0.base.ctype="S";
                v0.ctype="S";
                v=cellfun(@(p)v0.set_cp(1,p),mat2cell(eye(N),ones(1,N),N));
            end
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
                % assert(b==Uqsl2.B)
                ret=arr{1}(p);
            end
            function ret=createDelta(N)
                [O,E,F,K,Ki]=obj.getGenerator(N);
                I_=O.unit;
                % arr=[arrayfun(@(n)(E(n)|I_)+(K(n)|E(n)),1:N),...
                %     arrayfun(@(n)(F(n)|Ki(n))+(I_|F(n)),1:N),...
                %     arrayfun(@(n)K(n)|K(n),1:N),...
                %     arrayfun(@(n)Ki(n)|Ki(n),1:N),...
                %     I_];
                  arr=[arrayfun(@(n)(E(n)|K(n))+(I_|E(n)),1:N),...
                    arrayfun(@(n)(F(n)|I_)+(Ki(n)|F(n)),1:N),...
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