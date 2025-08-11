classdef(InferiorClasses=?sym) Usl2<strAlg&UEAlg
    properties(Constant,Hidden)
        B=Bases(3,["E" "F" "H"],"sl2")
    end
    
    %% generation
    methods(Static)
        function obj=make(cf,pw,~)
            obj=Usl2().make@strAlg(cf,pw);
        end
        function [O,E,F,H]=getGenerator(obj)
            O=Usl2();
            O=O.make(0,{[]});
            E=O.make(1,{1});
            F=O.make(1,{2});
            H=O.make(1,{3});
        end
    end
    methods
        
        %% relation
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                S.rel(1)=Usl2.make([1 -1 -1],{[1 2] [2 1] 3});
                S.rel(2)=Usl2.make([1 -1 -2],{[3 1] [1 3] 1});
                S.rel(3)=Usl2.make([1 -1 2],{[3 2] [2 3] 2});
                S.comm=[nan;nan];
                S.inv=[nan;nan];
                S=obj.get2vRelation_(S);
                % disp("更新")
            end
            rel=S.rel;
            mlist=S.mlist;
            comm=S.comm;
            inv=S.inv;
        end

        
    end
end