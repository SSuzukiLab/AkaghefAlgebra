classdef(InferiorClasses=?sym) StrWeylAlg<StrAlg&UEAlg
    properties(Constant,Hidden)
        B=TypeParam(@(N)Bases(2*N,["X"+(1:N) "D"+(1:N)],"weyl_"+N))
    end
    properties
        Nvar
       
    end
    
    %% generation
    methods
        function obj=StrWeylAlg(Nvar)
            obj.Nvar=Nvar;
            obj.algbase=obj.B.get(Nvar);
        end
        function obj=make(obj,cf,pw)
            obj=obj.make@StrAlg(cf,pw,obj.algbase);
        end
    end
    methods(Static)
        function [O,varargout]=getGenerator(Nvar,isNumberedVarName)
            arguments
                Nvar 
                isNumberedVarName =false;
            end
            O=StrWeylAlg(Nvar);
            O=O.make(0,{[]});
            C=num2cell(1:Nvar);
            C(2,:)=cellfun(@(i){O.make(1,{i})},C(1,:));
            C(3,:)=cellfun(@(i){O.make(1,{Nvar+i})},C(1,:));
            X=dictionary(C{[1 2],:});
            D=dictionary(C{[1 3],:});
            if isNumberedVarName
            varargout={X,D};
            else
                for i=1:Nvar
                    assignin('base',"X"+i, X(i));
                    assignin('base',"D"+i, D(i));
                end
                varargout={};
            end

        end
    end
    methods
        
        %% relation
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent RS
            if isempty(RS)
                RS=TypeParam(@createRel);
            end
            S=RS.get(obj.Nvar);
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
            function S_=createRel(Nvar)
                O=StrWeylAlg.getGenerator(Nvar);
                S_=Struct;
                % DX-XD=1
                S_.rel=arrayfun(@(i)O.make([1 -1 -1],{[Nvar+i i] [i Nvar+i] []}),1:Nvar);
                S_.comm=[nan;nan];
                S_.inv=[nan;nan];
                S_=O.get2vRelation_(S_);
            end
        end

        %% Hopf
        % comult
        function ret = Delta(obj)
            I=obj.unit;
            ret=obj.algfun(@fun,I|I);
            function ret=fun(p,b)
                % primitive-like
                X=I.set_cp(1,{p},{b});
                ret=(X|I)+(I|X);
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
            ret=obj.lfun(@fun);
            function [c,p,b]=fun(p,b)
                c=(-1)^length(p{1});
            end
        end

    end
end