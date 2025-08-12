classdef Usl2Borel<StrAlg&UEAlg
    properties(Constant)
        PairClass='FSL2Borel'
        B=Bases(2,["E" "H"],"sl2Borel")
    end
    methods(Static)
        function [O,E,H]=getGenerator
            O=Usl2Borel();
            O=O.make(0,{[]});
            E=O.make(1,{1});
            H=O.make(1,{2});
        end
    end
    methods
        function obj=make(obj,cf,pw,~)
            obj=obj.make@StrAlg(cf,pw,Usl2Borel.B);           
        end
        function ret=getP(obj,p)
            
        end

        %% relation
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                S.rel(1)=obj.make([1 -1 -2],{[2 1] [1 2] 1});
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
        function ret=HPGen(i1,i2)
            ret=HPGen(i2,i1);
        end
    end

end