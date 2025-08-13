classdef(InferiorClasses=?sym) StrC2modalg<StrAlg
    properties(Constant,Hidden)
        B=Bases(4,["x_m2" "x_m1" "x_1" "x_2"],"StrC2modalg")
    end

    methods
        
        function obj=StrC2modalg(varargin)
            obj@StrAlg(varargin{:});
            obj.priority=1:4;
            obj.ctype="S";
        end
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                J=dictionary([-2 -1 1 2],1:4);
                q=sym('q');
                % KE=q^2EK
                O=StrC2modalg.getGenerator;
                S.rel=O.make([1 -q],{J([-1 -2]) J([-2 -1])});
                S.rel(end+1)=O.make([-1 q],{J([1 -2]) J([-2 1])});
                S.rel(end+1)=O.make([-1 q],{J([2 1]) J([1 2])});
                S.rel(end+1)=O.make([-1 q],{J([2 -1]) J([-1 2])});
                S.rel(end+1)=O.make([-1 q^2],{J([2 -2]) J([-2 2])});
                S.rel(end+1)=O.make([-1 q^2 q^2*(q-q^-1)],{J([1 -1]) J([-1 1]) J([-2 2])});
                S.comm=[nan;nan];
                S.inv=[nan;nan];
                S=obj.get2vRelation_(S);
            end
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
        end
        function arg=casttype(obj,arg)
            arg=StrC2modalg(arg);
        end
        function obj=make(obj,cf,pw,~)
            obj=obj.make@StrAlg(cf,pw,StrC2modalg.B);
        end
    end
    methods(Static)
        
        function [O,Xm2,Xm1,X1,X2]=getGenerator()
            O=StrC2modalg().make(0,{});
            Var=["Xm2","Xm1","X1","X2"];
            for i=1:4
                X=O.make(1,{i});
                eval(Var(i)+"=X;");
            end
        end
    end

end