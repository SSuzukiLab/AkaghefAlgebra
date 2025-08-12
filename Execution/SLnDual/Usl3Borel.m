classdef Usl3Borel<StrAlg&UEAlg
    properties(Constant)
        PairClass='FSL3Borel'
        B=Bases(2,["E1", "E12", "E2", "H1", "H2"],"sl3Borel")
    end
    methods(Static)
        function [O,E,H]=getGenerator
            O=Usl3Borel();
            O=O.make(0,{[]});
            E=dictionary( ...
                1,O.make(1,{1}), ...
                12,O.make(1,{2}), ...
                2,O.make(1,{3}));
            H=dictionary(...
                1,O.make(1,{4}), ...
                2,O.make(1,{5}));
        end
    end

    methods

        function obj=make(obj,cf,pw,~)
            obj=obj.make@StrAlg(cf,pw,Usl3Borel.B);
            
        end
        function obj=getP(obj,p)
            obj.make(1,{p});
        end

        %% relation
        function [rel,mlist,comm,inv]=get2vRelation(obj)
            persistent S
            if isempty(S)
                S=struct;
                % S.rel(1)=obj.make([1 -1 -2],{[2 1] [1 2] 1});
                S.rel=obj.empty;
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